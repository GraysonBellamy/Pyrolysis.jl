# 16. Conservation Diagnostics

This chapter derives the runtime conservation diagnostics that Pyrolysis.jl uses
to audit its own numerics and to quantify the physical cost of the Formulation A
modeling approximations. We develop a per-component **mass ledger**, a **dual
energy ledger** (a *discrete* form that closes to integrator tolerance by
construction and a *physical* form that exposes the cost of the quasi-steady-gas
approximation), and a **per-reaction mass-flow ledger**. For every accumulator we
give the continuous balance law, its exact discrete counterpart as the code
assembles it, the precise sign conventions, and the conditions under which the
residual is expected to vanish. The chapter closes by stating exactly which
balance terms are *not* yet accumulated, why the discrete-form residual cannot
close for time-varying composition without them, and how the production
diagnostics callback operates in snapshot mode.

All symbols follow the unified nomenclature of §2 (Nomenclature). In this chapter
`i` is the **cell index** and `r` the **reaction index**; `j` is the **component
index** (1…$N_C$). The three binding sign conventions are in force throughout:
flux is positive in the $+z$ direction (bottom→top); the discrete divergence
$(F_R-F_L)A/V$ is positive for *outflow*; and the heat-of-reaction convention is
$h>0$ endothermic, $h<0$ exothermic, with the volumetric heat source
$Q_\mathrm{rxn}=-\sum_r h_r r_r$ (see §6, Reaction Kinetics, and §7, Heat
Conduction).

The implementing module is `Diagnostics`
(`src/Diagnostics/Diagnostics.jl`), comprising the conservation-check kernels
(`src/Diagnostics/conservation_checks.jl`) and the callback infrastructure
(`src/Diagnostics/callback.jl`).

---

## 16.1 Purpose and design philosophy

A pyrolysis solver couples a stiff reaction network, in-depth radiation, gas
transport, a moving (ALE) mesh, and adaptive refinement, all advanced by an
implicit time integrator. With this many interacting mechanisms, a silent error
in flux or source assembly can masquerade as physics. The diagnostics subsystem
exists to make conservation **observable at runtime** rather than only in offline
post-processing, so that a broken term announces itself as a non-closing ledger.

Two distinct questions must be separated, and the design answers them with **two
different energy ledgers**:

1. *Is the discretized RHS being integrated correctly?* This is a question of
   bookkeeping: do the accumulated boundary, reaction, and gas-coupling terms
   add up to the change in stored energy that the RHS actually evolves? The
   **discrete-form ledger** answers this and is designed to close to integrator
   tolerance *by construction* — it is assembled from exactly the terms the RHS
   computes. A failure here is a bug in flux/source assembly or in the callback
   plumbing, not a physics error.

2. *What does the Formulation A approximation cost?* Formulation A (§3, §7)
   removes gas sensible-heat storage from the volumetric heat capacity
   $\rho c_p^\mathrm{eff}$ and carries gas enthalpy entirely through the
   advection source $S_\mathrm{conv}$ (the quasi-steady-gas approximation). The
   **physical-form ledger** includes the *total* sensible storage
   $\rho c_p^\mathrm{total}$ (matrix + gas) and quantifies, as a residual, how
   much energy the matrix equation fails to account for. For canonical pyrolysis
   this residual is below 1 % of the initial energy; if it grows large, the
   regime is stressing the approximation and the user is warned.

The key design subtlety, stated plainly in the source documentation
(`conservation_checks.jl:243–247`), is that the discrete-form ledger audits
*bookkeeping, not the physical model*: because it is built from the terms the RHS
actually assembles, **it will close even if the RHS itself double-counts a
mechanism.** That is precisely why the second, independent physical-form ledger
exists, and why the per-reaction mass ledger is computed independently of the
global mass ledger.

> **Implementation note.** The three trackers and the callback live in module
> `Diagnostics`. They are pure functions of the *current mesh state*
> (`mesh.cell_states`, `mesh.cell_volumes`, `mesh.active_cell_mask`) plus a small
> set of mutable accumulators; they never read the raw state vector directly. The
> solver synchronizes the mesh from the integrator's `u` before each diagnostic
> snapshot via `update_mesh_from_state!`
> (`src/Solver/solve.jl:886`).

---

## 16.2 The continuous balance laws being audited

Before discretizing, we state the integral conservation laws that the ledgers
approximate. Let $\Omega(t)=[0,L(t)]\times A(t)$ be the (possibly moving,
variable-area) domain.

### 16.2.1 Mass

Integrating the component conservation law (§3.2) over the domain and applying
the divergence theorem, for each component $j$:

$$
\frac{d}{dt}\int_\Omega \xi_j \,dV
= -\int_{\partial\Omega} J_j\cdot \hat{\mathbf n}\,dA
\;+\; \int_\Omega \Big(\tfrac{\partial \xi_j}{\partial t}\Big)_\mathrm{rxn}\,dV .
$$

The reaction term integrates to zero *summed over all components* because the
single-reactant stoichiometry conserves mass exactly ($\sum_j \nu_{r,j}=0$, with
the reactant coefficient $-1$ balancing product yields summing to $1$; see §6.5).
Therefore the **total** mass changes only through the boundary efflux. For a
component-by-component statement we keep the reaction term explicit but rely on
build-time stoichiometric closure (§6.8) to guarantee $\sum_j S_j=0$. With the
two boundaries (top $z=L$, bottom $z=0$) and the sign convention that $J_j>0$ is
$+z$ transport, the *outward* normal flux is $+J_j$ at the top and $-J_j$ at the
bottom, so the net mass leaving is $(J_{\mathrm{top},j}-(-J_{\mathrm{bot},j}))$
times the area. The diagnostics layer instead adopts the **into-domain**
convention for its boundary inputs (positive $J$ = inflow), under which the net
efflux is $-(J_{\mathrm{top},j}+J_{\mathrm{bot},j})A$ — see §16.3.

### 16.2.2 Energy

Multiplying the conservative-form energy equation by the control volume and
integrating, the *total sensible* energy obeys

$$
\frac{d}{dt}\int_\Omega \rho c_p^\mathrm{total}\,T\,dV
= \underbrace{\int_{\partial\Omega} q_\mathrm{BC}\,dA}_{\text{boundary input}}
+ \underbrace{\int_\Omega Q_\mathrm{rxn}\,dV}_{\text{reaction}}
- \underbrace{\sum_g \int_{\partial\Omega}\dot m_g\,\Delta h_g(T_s,T_\infty)\,dA}_{\text{enthalpy carried out by gas}} ,
$$

where $\rho c_p^\mathrm{total}=\sum_{j} \xi_j c_{p,j}(T)$ sums **all** components
(matrix + gas). This is the law the physical-form ledger discretizes.

Formulation A, however, evolves a *non-conservative* temperature equation with a
*matrix-only* storage term (§7.3),
$\rho c_p^\mathrm{eff}=\sum_{j\in\{s,\ell\}}\xi_j c_{p,j}(T)$, and carries gas
enthalpy through the volumetric sources $S_\mathrm{conv}$ (always on) and
$S_\mathrm{gen}$ (optional). The energy the *matrix equation* actually evolves
therefore obeys

$$
\frac{d}{dt}\int_\Omega \rho c_p^\mathrm{eff}\,T\,dV
\;\approx\; \int_{\partial\Omega} q_\mathrm{cond}\,dA
+ \int_\Omega Q_\mathrm{rxn}\,dV
+ \int_\Omega S_\mathrm{conv}\,dV
+ \int_\Omega S_\mathrm{gen}\,dV ,
$$

with $q_\mathrm{cond}$ the *conducted* boundary flux (not the full applied BC
flux). This is the law the discrete-form ledger discretizes. The $\approx$ flags
two terms missing from the storage time-derivative — the composition-change term
and the ALE volume-change term — which are the subject of §16.9.

The difference between the two storage definitions,

$$
\Delta E_\mathrm{gas}(t) \equiv E_\mathrm{total}(t) - E_\mathrm{matrix}(t)
= \int_\Omega \Big(\sum_{g}\xi_g c_{p,g}(T)\Big) T\,dV ,
$$

is the **gas-storage gap** — the sensible heat physically held in the in-pore gas
that the matrix equation never evolves. Its smallness ($\ll 1\%$ at 1 atm) is the
justification for Formulation A, and the physical-form residual is the runtime
audit of that smallness.

---

## 16.3 Mass conservation ledger

### 16.3.1 Domain mass

The total mass of component $j$ in the domain at the current state is the
volume-weighted sum of cell concentrations over the **active** cells only:

$$
m_j(t) = \sum_{i\in\mathrm{active}} \xi_{j,i}\,V_i .
$$

This is a direct first-order (midpoint) quadrature of $\int_\Omega \xi_j\,dV$:
$\xi_{j,i}$ is the cell-centered mass concentration [kg·m⁻³] and $V_i$ the cell
volume [m³] (`mesh.cell_volumes[i]`, which is time-dependent in ALE mode and
carries the variable cross-section). The active-cell mask
(`mesh.active_cell_mask[i]`) excludes cells that AMR has refined away (see §14,
Adaptive Mesh Refinement), so the ledger always reflects only the live solution.

> **Implementation note.** `compute_domain_mass` /
> `src/Diagnostics/conservation_checks.jl:54–71`. The loop iterates
> `eachindex(mesh.cells)` and `continue`s on inactive cells; it accumulates into
> a `MVector{NC,Float64}` and returns an `NTuple{NC,Float64}`. The mask test must
> be performed cell-by-cell — never `eachindex` over a compacted array — which is
> the AMR invariant noted in §14.6.

### 16.3.2 Boundary efflux accumulation

The boundary mass flux update integrates the surface fluxes over one substep.
With the diagnostics **into-domain** sign convention ($J>0$ = inflow), the mass
that *leaves* during a substep of size $\Delta t$ is

$$
\Delta m_{\mathrm{out},j} = -\big(J_{\mathrm{top},j} + J_{\mathrm{bot},j}\big)\,A\,\Delta t ,
$$

accumulated into a running `cumulative_efflux`:

$$
m_{\mathrm{efflux},j}(t) = \int_0^t -\big(J_{\mathrm{top},j}+J_{\mathrm{bot},j}\big)\,A\,d\tau
\;\;\approx\;\; \sum_{\text{substeps}} \Delta m_{\mathrm{out},j}.
$$

A positive flux (inflow) contributes *negatively* to the efflux (mass stays in
the domain); a negative flux (the typical case — gas escaping at the top,
$J_{\mathrm{top},j}<0$) contributes *positively*. The accumulator therefore
counts up the total mass that has left. The canonical configuration is
$J_{\mathrm{top},j}<0$ (gas escaping) and $J_{\mathrm{bot},j}=0$ (impermeable
back face), giving a monotonically increasing efflux for the gaseous products.

> **Implementation note.** `update_mass_tracker!` /
> `conservation_checks.jl:94–116`. `J_boundary` is an `NC×2` matrix with column 1
> = top (boundary id 1), column 2 = bottom (boundary id 2); the code computes
> `mass_out = -(J_top + J_bottom) * A * dt` per component and adds it to
> `cumulative_efflux`. It also advances `tracker.last_time += dt`.

**Variable cross-section (ALE).** For a moving mesh the area $A(t)$ changes within
the substep. To keep the boundary contribution consistent with the trapezoidal
geometry update used elsewhere (§11, ALE Framework), the caller passes the
**trapezoidal-in-area average** over the substep,

$$
A = \tfrac12\big(A_t^\mathrm{pre} + A_t^\mathrm{post}\big),
$$

so that $\Delta m_{\mathrm{out},j}$ uses the substep-mean area. With the identity
lateral-shrinkage law (`NoChi`) or a non-ALE (Eulerian) run, $A$ is constant and
equals `mesh.cross_section_area_initial` (`conservation_checks.jl:83–87`).

### 16.3.3 The balance check and tolerance

Mass conservation is the statement that the current mass equals the initial mass
minus everything that has escaped:

$$
m_j(t) \;\overset{?}{=}\; m_{\mathrm{initial},j} - m_{\mathrm{efflux},j}(t)
\quad\Longleftrightarrow\quad
m_j(t) + m_{\mathrm{efflux},j}(t) = m_{\mathrm{initial},j}.
$$

The check defines the **expected mass** $m_{\mathrm{expected},j}=m_{\mathrm{initial},j}-m_{\mathrm{efflux},j}$
and a per-component relative error

$$
\epsilon_{m,j} =
\begin{cases}
\dfrac{\big|m_j(t) - m_{\mathrm{expected},j}\big|}{\,|m_{\mathrm{initial},j}|\,},
& |m_{\mathrm{initial},j}| > \texttt{eps()},\\[2.2ex]
\big|m_j(t) - m_{\mathrm{expected},j}\big|, & \text{otherwise (absolute fallback).}
\end{cases}
$$

The **zero-initial-mass fallback** is essential: gaseous products typically start
with $m_{\mathrm{initial},j}\approx 0$, so a relative error would divide by zero.
Switching to absolute error correctly flags even a tiny *spurious* mass when
there was no baseline to normalize against. The verdict is

$$
\text{conserved} \iff \max_j \epsilon_{m,j} < \texttt{tol},\qquad \texttt{tol}=10^{-10}.
$$

This is a tight, *bookkeeping-level* tolerance: it should hold to roundoff
because both sides are computed by exactly consistent quadratures. The check
returns the full per-component breakdown (`current_mass`, `expected_mass`,
`relative_error`, `max_error`, `conserved`).

**Scope of the $10^{-10}$ tolerance.** The roundoff-level closure holds where
the discrete ledger telescopes exactly *and* every escape channel is counted:
condensed components (mesh-bound, flux-form transport), all components on a
fixed (Eulerian) mesh, and gas components on a **fully closed** ALE column —
the gas mesh-motion correction is the conservative flux form (§13.8.3), so a
collapsing impermeable column conserves each gas component to integrator
tolerance (V-CPL-03 P4). The remaining gap is bookkeeping, not physics: with
an *open* mass BC on a moving boundary, the boundary face additionally sweeps
gas in/out at $\xi_{\text{face}}\,w_{\text{face}}A$ (a regressing surface
releases the pore gas it uncovers), and the tracker's expected-mass formula
integrates only the transport flux $J_{\text{gas}}$ — expect a corresponding
offset in open-boundary ALE runs.

> **Implementation note.** `check_mass_conservation` /
> `conservation_checks.jl:145–177`; `mass_conservation_report`
> (lines 187–218) renders the per-component table and an overall
> PASSED/FAILED line. The state is held in `MassConservationTracker{NC}`
> (`conservation_checks.jl:28–32`): `initial_mass`,
> `cumulative_efflux`, `last_time`.

### 16.3.4 Sources of (small) mass-ledger discrepancy

Even when the implementation is correct the ledger carries two avoidable and one
unavoidable error source. (i) **Offline vs. online integration.** If the efflux
is reconstructed offline by trapezoidal integration of saved surface fluxes
(rather than accumulated online per substep), the integration error of the
post-processing quadrature appears in the residual — the `energy_ledger_audit.jl`
example explicitly reconstructs $\int \dot m_g\,A\,dt$ this way and expects
$\ll 0.1\%$ (`examples/inverse_analysis/energy_ledger_audit.jl:174–184`).
(ii) **Substep granularity.** The online accumulator uses the substep-mean area;
on a rapidly shrinking ALE mesh the trapezoidal-area approximation is
second-order in $\Delta t$. (iii) **Active-mask transitions.** When a surface
cell is merged or refined away (§14.3), its mass is redistributed by an
energy/mass-conserving remap; the ledger picks up the post-remap masses, so a
correctly conserving remap leaves no residual but a buggy one would surface here.
(iv) **ALE gas-species advection (unavoidable at fixed resolution).** The
$+w\,\partial_z\xi_g$ mesh-motion correction is applied pointwise
(non-conservative form), so under ALE the gas-species ledger residual is
truncation-order rather than roundoff — see the tolerance-scope note in
§16.3.3. Solids close exactly regardless.

---

## 16.4 Energy storage: the two heat-capacity definitions

The two energy ledgers differ only in which volumetric heat capacity weights the
stored energy. Both are computed by the same `@generated` summation kernel
(`effective_properties.jl:305–320`), branching on the compile-time
`is_gas` predicate so the branch is eliminated at codegen.

**Matrix-only (gas-excluded) energy** — the storage the matrix equation evolves:

$$
E_\mathrm{matrix}(t) = \sum_{i\in\mathrm{active}} \rho c_p^\mathrm{eff}(\xi_i, T_i)\, T_i\, V_i,
\qquad
\rho c_p^\mathrm{eff}(\xi,T)=\sum_{j\in\{\text{solid},\text{liquid}\}}\xi_j\,c_{p,j}(T).
$$

**Total energy** — matrix plus all gas storage:

$$
E_\mathrm{total}(t) = \sum_{i\in\mathrm{active}} \rho c_p^\mathrm{total}(\xi_i, T_i)\, T_i\, V_i,
\qquad
\rho c_p^\mathrm{total}(\xi,T)=\sum_{j=1}^{N_C}\xi_j\,c_{p,j}(T).
$$

Both sums respect the active-cell mask. The product $\rho c_p\cdot T$ uses the
**absolute** temperature $T$ (not $T-T_\mathrm{ref}$): the diagnostics measure a
sensible-energy *content* relative to $0\,\mathrm{K}$, which is acceptable because
the ledgers only ever use *differences* $E(t)-E(0)$ in which the datum cancels
(provided $\rho c_p$ is held fixed; see the composition-change caveat in §16.9).

> **Implementation note.** `compute_sensible_energy` /
> `conservation_checks.jl:306–322` (matrix) and
> `compute_total_sensible_energy` / lines 333–349 (total). They call
> `Properties.effective_heat_capacity` and `Properties.total_heat_capacity`
> respectively. The two heat-capacity functions are
> `effective_properties.jl:341–343` and `356–369`; note the *only*
> structural difference is the `if !is_gas(...)` guard in the matrix version.

### 16.4.1 The gas-storage gap

The difference

$$
\Delta E_\mathrm{gas}(t) = E_\mathrm{total}(t) - E_\mathrm{matrix}(t)
= \sum_{i\in\mathrm{active}}\Big(\sum_{g}\xi_{g,i}\,c_{p,g}(T_i)\Big)\,T_i\,V_i \;\ge 0
$$

is the sensible energy held in the in-pore gas. It is the *physical* magnitude of
the quasi-steady-gas approximation. Order-of-magnitude: gas concentrations in the
pore space are $O(\rho_g)\sim P M/(R_g T)$, which at 1 atm and 700 K with
$M\sim0.03$ kg·mol⁻¹ is $\sim0.5$ kg·m⁻³, against matrix densities of hundreds of
kg·m⁻³ — hence $\Delta E_\mathrm{gas}/E_\mathrm{matrix}\sim10^{-3}$. The gap
scales with pressure (more confined gas mass), so it is the diagnostic to watch
for high-pressure or Darcy-dominated runs. The runtime row reports both
$\Delta E_\mathrm{gas}$ in joules and, in the audit example, its fraction of
$E_\mathrm{matrix}$ (`energy_ledger_audit.jl:170–172`).

---

## 16.5 Discrete-form energy ledger

### 16.5.1 Continuous-to-discrete construction

The discrete-form ledger is assembled from *exactly the terms the RHS computes*.
Integrating the matrix-energy balance of §16.2.2 in time and discretizing each
term by the same quadrature the RHS uses gives

$$
E_\mathrm{matrix}(t) = E_\mathrm{matrix}(0)
+ \int_0^t \Big[\, q_\mathrm{cond,boundary}\,A
+ Q_\mathrm{rxn}^\mathrm{tot}
+ S_\mathrm{conv}^\mathrm{tot}
+ S_\mathrm{gen}^\mathrm{tot} \,\Big]\,d\tau,
$$

where each integrand is the *volume-integral of the corresponding RHS field*:

- $q_\mathrm{cond,boundary}$ — net **conducted** boundary heat flux [W·m⁻²],
  positive into the slab. This is the flux the Newton boundary solver actually
  drives into the first/last cell (§7.6, §12.2), **not** the externally applied
  BC flux.
- $Q_\mathrm{rxn}^\mathrm{tot}=\sum_i Q_{\mathrm{rxn},i}V_i = \int_\Omega Q_\mathrm{rxn}\,dV$ [W],
  the total reaction heat ($Q_\mathrm{rxn}=-\sum_r h_r r_r$; §6.6).
- $S_\mathrm{conv}^\mathrm{tot}=\int_\Omega S_\mathrm{conv}\,dV$ [W], the total gas
  advective enthalpy sink; **typically negative** (heat carried out by gas; §7.4).
- $S_\mathrm{gen}^\mathrm{tot}=\int_\Omega S_\mathrm{gen}\,dV$ [W], the optional
  generation sink; identically zero unless `energy_form = :with_generation_sink`
  (§7.5).

### 16.5.2 Per-substep accumulation

The energy tracker accumulates the time-integrals as Euler products over each
substep $\Delta t$:

$$
\begin{aligned}
\texttt{cum\_boundary\_conduction} &\mathrel{+}= q_\mathrm{cond,boundary}\,A\,\Delta t,\\
\texttt{cum\_boundary\_applied}    &\mathrel{+}= q_\mathrm{BC,applied}\,A\,\Delta t,\\
\texttt{cum\_reaction\_heat}       &\mathrel{+}= Q_\mathrm{rxn}^\mathrm{tot}\,\Delta t,\\
\texttt{cum\_advection\_sink}      &\mathrel{+}= S_\mathrm{conv}^\mathrm{tot}\,\Delta t,\\
\texttt{cum\_generation\_sink}     &\mathrel{+}= S_\mathrm{gen}^\mathrm{tot}\,\Delta t,\\
\texttt{cum\_transpiration}        &\mathrel{+}= q_\mathrm{transp}\,A\,\Delta t .
\end{aligned}
$$

The boundary terms carry the area factor $A$ and have flux units [W·m⁻²], so
$q\cdot A\cdot\Delta t$ has units of joules; the source terms are already
volume-integrals [W], so $\cdot\,\Delta t$ gives joules. Both
`cum_boundary_conduction` and `cum_boundary_applied` are accumulated every
substep so the tracker can serve both ledgers; the discrete ledger uses
*conduction* and the physical ledger uses *applied* (§16.6).

> **Implementation note.** `update_energy_tracker!` /
> `conservation_checks.jl:372–391`. The accumulators live in
> `EnergyConservationTracker` (lines 275–285): `initial_energy_matrix`,
> `initial_energy_total`, `cumulative_boundary_conduction`,
> `cumulative_boundary_applied`, `cumulative_reaction_heat`,
> `cumulative_advection_sink`, `cumulative_generation_sink`,
> `cumulative_transpiration`, `last_time`.

### 16.5.3 The discrete residual

At a check point the *expected* matrix energy is the initial energy plus the
accumulated discrete sources (note: **conduction**, **not** applied, and **with**
both $S_\mathrm{conv}$ and $S_\mathrm{gen}$):

$$
E_\mathrm{matrix}^\mathrm{exp}(t) = E_\mathrm{matrix}(0)
+ \texttt{cum\_boundary\_conduction}
+ \texttt{cum\_reaction\_heat}
+ \texttt{cum\_advection\_sink}
+ \texttt{cum\_generation\_sink},
$$

and the **discrete residual** is the gap between the measured and expected stored
energy:

$$
\boxed{\,r_\mathrm{discrete}(t) = E_\mathrm{matrix}(t) - E_\mathrm{matrix}^\mathrm{exp}(t)\,}
$$

with relative magnitude $r_\mathrm{discrete}/|E_\mathrm{matrix}(0)|$ (absolute if
$|E_\mathrm{matrix}(0)|\le\texttt{eps()}$) and the verdict
$r_\mathrm{discrete}^\mathrm{rel}<\texttt{tol}$, default $\texttt{tol}=10^{-6}$.

> **Implementation note.** `check_energy_conservation` /
> `conservation_checks.jl:412–463`. The expected-matrix expression is lines
> 423–428 and `discrete_residual = E_matrix_now - expected_matrix` on line 429.
> The relative form (lines 440–444) divides by `|initial_energy_matrix|` when
> non-trivial.

### 16.5.4 Why it closes — and what it does *not* prove

This residual closes to integrator tolerance *by construction* whenever the
heat capacity and volume are constant over the trajectory, because both
$E_\mathrm{matrix}(t)$ and its expectation are assembled from the same RHS fields
under exactly consistent quadrature; the only error is the integrator's local
truncation error in advancing $T$, which the implicit solver controls to
`reltol`. A non-closing discrete residual therefore points to a **bookkeeping
bug**: a flux computed but not accumulated, a callback firing at the wrong time,
a sign error in a source term, or a missing factor of $A$ or $V$.

The crucial limitation, repeated from §16.1: the discrete ledger **cannot detect
a double-count**. If the RHS were (incorrectly) to add both $S_\mathrm{conv}$ and
a surface transpiration enthalpy term, the accumulator would faithfully record
both, the expectation would include both, and the residual would still close.
Catching that class of error is the job of the physical-form ledger (§16.6) and
of the construction-time mutual-exclusivity guard between $S_\mathrm{conv}$ and
the transpiration BC (§7.6). It is also subject to the composition-change and
ALE volume-change omissions of §16.9, which prevent *exact* closure whenever the
mixture heat capacity or cell volume vary in time.

---

## 16.6 Physical-form energy ledger

### 16.6.1 The total-energy balance

The physical-form ledger discretizes the *conservative total-energy* law of
§16.2.2. The expected total energy is

$$
E_\mathrm{total}^\mathrm{exp}(t) = E_\mathrm{total}(0)
+ \texttt{cum\_boundary\_applied}
+ \texttt{cum\_reaction\_heat}
- \texttt{cum\_transpiration},
$$

and the **physical residual** is

$$
\boxed{\,r_\mathrm{physical}(t) = E_\mathrm{total}(t) - E_\mathrm{total}^\mathrm{exp}(t)\,}
$$

with relative magnitude $r_\mathrm{physical}/|E_\mathrm{total}(0)|$ and verdict
$r_\mathrm{physical}^\mathrm{rel}<\texttt{tol}$.

Three deliberate differences from the discrete ledger encode the physics:

1. **Total storage.** It uses $E_\mathrm{total}$ (with gas), so the gas-storage
   gap appears as a genuine, physically meaningful contribution rather than
   being silently excluded.
2. **Applied (not conducted) boundary flux.** It uses `cum_boundary_applied`,
   built from $q_\mathrm{BC,applied}$ — the *externally applied* driving flux
   (which may include incident radiation and convective exchange), not the
   conducted flux that the matrix equation sees. This represents the true energy
   crossing the boundary.
3. **Explicit transpiration sink.** It subtracts the cumulative transpiration —
   the sensible enthalpy physically carried out of the domain by escaping gas at
   the surface temperature — as energy *leaving*.

### 16.6.2 The transpiration sign

The transpiration term enters $E_\mathrm{total}^\mathrm{exp}$ with a **minus**
sign, equivalently it enters $r_\mathrm{physical}$ with a **plus** sign. This is
counterintuitive but correct. The cumulative transpiration

$$
\texttt{cum\_transpiration}(t) = \int_0^t \sum_g \dot m_g(T_s)\,\Delta h_g(T_s,T_\infty)\,A\,d\tau \;\ge 0
$$

is the magnitude of sensible enthalpy that escaping gas carries away. Energy that
*leaves* must be subtracted from the energy *delivered* to predict the remaining
stored energy:

$$
E_\mathrm{total}(t) \;\overset{?}{=}\; \underbrace{E_\mathrm{total}(0)
+ \texttt{cum\_boundary\_applied} + \texttt{cum\_reaction\_heat}}_{\text{energy in}}
\;-\; \underbrace{\texttt{cum\_transpiration}}_{\text{energy out via gas}} .
$$

Rearranging to the residual form $r_\mathrm{physical}=E_\mathrm{total}(t)-E_\mathrm{total}^\mathrm{exp}(t)$
puts $+\,\texttt{cum\_transpiration}$ in the residual (`conservation_checks.jl:434–438`,
and the rationale at lines 476–484 of the source).

### 16.6.3 Conservative-form recovery and telescoping

When does the physical residual vanish exactly? When the *full conservative form*
is recovered. The gas-coupled volumetric divergence splits as (§7.4–7.5)

$$
\nabla\!\cdot\!\Big(\sum_g c_{p,g}(T-T_\mathrm{ref})\,J_g\Big)
= \underbrace{\sum_g c_{p,g}\,J_g\!\cdot\!\nabla T}_{-S_\mathrm{conv}}
+ \underbrace{(T-T_\mathrm{ref})\sum_g c_{p,g}\,\nabla\!\cdot\!J_g}_{-S_\mathrm{gen}\;(\text{for const }c_p)} ,
$$

and for the *true* sensible enthalpy $\Delta h_g$ the divergence-theorem identity
over the domain reads

$$
\int_\Omega \big(S_\mathrm{conv}+S_\mathrm{gen}\big)\,dV
= -\sum_g \dot m_{g,\mathrm{surface}}\,\Delta h_g(T_s,T_\infty).
$$

This is the **telescoping identity**: the volumetric advection plus generation
sink, integrated over the domain, exactly equals the *negative* of the surface
transpiration outflow. It holds *exactly* in the discrete scheme because
$S_\mathrm{gen}$ uses the same midpoint-rule quadrature for $\Delta h_g$ and the
same reference temperature $T_\mathrm{ref}=T_\infty$ as the transpiration BC (§7.5).

The consequences:

- In **`:with_generation_sink`** mode *paired with* the surface transpiration BC,
  the conservative form is recovered to FVM truncation error: the physical
  residual closes (modulo the gas-storage gap evolution and the §16.9
  omissions).
- In the default **`:standard`** mode, $S_\mathrm{gen}$ is dropped, so only
  $\int S_\mathrm{conv}$ telescopes; the volumetric integration of the
  gas-enthalpy divergence is closed *only at the boundary*, not in-depth. The
  physical residual is then expected to be *small but nonzero* — its magnitude is
  the deficit $\int S_\mathrm{gen}$ plus the gas-storage gap evolution, which is
  exactly the "cost" the ledger is built to expose.
- The transpiration BC must **never** be combined with $S_\mathrm{conv}$ alone:
  doing so applies the full surface outflow $G$ on top of the already-present
  $\int S_\mathrm{conv}$, double-counting the advected enthalpy by
  $\int S_\mathrm{conv}$. This is rejected at problem construction
  (`Problem.to_problem_def`; see §7.6 and the test
  `test/discretization/test_transpiration_energy_balance.jl`).

### 16.6.4 Expected magnitude and interpretation

For canonical pyrolysis (gas residence time $\ll$ thermal diffusion time, i.e.
1D slabs near atmospheric pressure), the physical residual should remain **below
1 % of the initial total energy**. A larger value means the regime is stressing
the quasi-steady-gas approximation — for example a high-pressure, low-porosity,
or strongly Darcy-driven case where in-pore gas storage and its transport timing
matter — and the source documentation directs the user to reconsider the setup
or enable `:with_generation_sink` with the transpiration BC
(`conservation_checks.jl:249–262`).

> **Implementation note.** The expected-total expression and physical residual
> are `conservation_checks.jl:434–438`; the relative form (lines 446–450)
> normalizes by `|initial_energy_total|`. The human-readable
> `energy_conservation_report` (lines 473–574) prints both ledgers side by side,
> including the gas-storage gap, with ✓/✗ verdicts.

---

## 16.7 Boundary-flux semantics: conducted vs. applied

The two energy ledgers depend on a careful distinction between two boundary
heat-flux quantities that coincide only for the simplest BCs.

- **Conducted flux $q_\mathrm{cond,boundary}$.** The heat actually conducted into
  the first/last interior cell, i.e. $k_\mathrm{eff}(T_s-T_\mathrm{cell})/\Delta z_{1/2}$
  evaluated at the converged surface temperature $T_s$ from the Newton boundary
  solve (§12.2). This is what the *matrix* energy equation receives, so the
  **discrete ledger** uses it.
- **Applied flux $q_\mathrm{BC,applied}$.** The externally applied driving flux
  that the boundary condition imposes — for a `RadiativeFluxBC` this includes the
  absorbed incident radiation and the re-radiation/convective exchange; for a
  pure `HeatFluxBC` it is the prescribed $q(t)$. This is the *true* energy
  crossing the surface, so the **physical ledger** uses it.

At steady state, or for an adiabatic/Dirichlet boundary, $q_\mathrm{cond}$ and
$q_\mathrm{BC,applied}$ coincide. They differ transiently whenever the surface
energy balance routes part of the applied flux into surface re-radiation,
convection, or transpiration rather than into conduction. The surface energy
balance the Newton solver enforces is (§12.2)

$$
q_\mathrm{BC,applied}(T_s) = q_\mathrm{cond}(T_s) + q_\mathrm{transp}(T_s),
$$

so that the difference $q_\mathrm{BC,applied}-q_\mathrm{cond}$ is precisely the
transpiration outflow at the surface (when the transpiration BC is active) or
the radiative/convective surface losses already folded into $q_\mathrm{BC,applied}$.
The `energy_ledger_audit.jl` example checks this identity directly as a
*per-time* surface-balance residual, expecting it to hold to the Newton tolerance
($\ll 10^{-7}$ relative) over the whole trajectory
(`energy_ledger_audit.jl:157–165, 192–203`). The solver exposes these as
`extras.surface_heat_flux_net` ($q_\mathrm{BC,applied}$ — the applied face
flux in every radiation regime; under `BEER_LAMBERT` it excludes the incident
radiation, which enters volumetrically),
`extras.surface_heat_flux_conduction` ($q_\mathrm{cond}$), and
`extras.surface_heat_flux_transpiration` ($q_\mathrm{transp}$), plus
`extras.volumetric_absorption` (the Beer–Lambert in-depth deposit
$\alpha_{\text{surf}}q_{\text{ext}}(1-e^{-\tau_{\text{tot}}})$; zero for other
models) and `extras.back_heat_flux` for the bottom boundary
(`src/Problem/pyrolysis_problem.jl`).

> **Implementation note.** Both fluxes feed `update_energy_tracker!` as separate
> keyword arguments (`q_cond_boundary`, `q_bc_applied`), are accumulated into
> separate cumulative fields, and are routed to the two ledgers respectively in
> `check_energy_conservation` (`conservation_checks.jl:423–438`).

---

## 16.8 Per-reaction mass-flow ledger

The global mass ledger (§16.3) verifies that the *domain* conserves mass. A
*finer* and *independent* check is that each reaction conserves mass on its own —
that products generated equal reactant consumed for every reaction. This catches
stoichiometry errors that the global ledger would miss (because mass consumed by
one reaction and generated by another can globally cancel).

### 16.8.1 Reaction extent and per-reaction accumulation

For reaction $r$ with cell-uniform volumetric rate $r_r$ [kg·m⁻³·s⁻¹], the mass
reacted over the active domain in one substep is

$$
m_{\mathrm{flow},r} = r_r\,V_\mathrm{total}\,\Delta t,
\qquad V_\mathrm{total}=\sum_{i\in\mathrm{active}} V_i .
$$

(The tracker takes a single representative or domain-averaged rate per reaction
and multiplies by the total active volume; for spatially uniform rates this is
exact, otherwise it is a domain-average estimate — see the open issue in
§16.8.3.) The accumulators update as

$$
\begin{aligned}
\texttt{reaction\_counts}_r &\mathrel{+}= m_{\mathrm{flow},r}
&&\text{(cumulative reaction extent }\textstyle\int r_r\,dt\,V_\mathrm{total}\text{),}\\
\texttt{reactant\_consumed}_{j_1(r),\,r} &\mathrel{+}= m_{\mathrm{flow},r}
&&\text{($j_1$ = first reactant of }r\text{; coefficient }1\text{),}\\
\texttt{product\_generated}_{j,\,r} &\mathrel{+}= \nu_{r,j}\,m_{\mathrm{flow},r}
&&\text{for each product }j\text{ of }r .
\end{aligned}
$$

The single-reactant convention (§6.5) makes the consumption side unambiguous: the
first reactant is consumed with implicit coefficient $1$, and each product is
generated with its explicit yield $\nu_{r,j}\ge0$.

> **Implementation note.** `update_reaction_tracker!` /
> `conservation_checks.jl:636–672`. It recomputes $V_\mathrm{total}$ over the
> active mask each call, then for each reaction adds `mass_flow` to the first
> reactant's `reactant_consumed[r1_idx, r]` and `ν * mass_flow` to each product's
> `product_generated[p_idx, r]`. State in `ReactionMassFlowTracker{NC,NR}`
> (lines 597–602): `reactant_consumed` [NC×NR], `product_generated` [NC×NR],
> `reaction_counts` [NR], `last_time`.

### 16.8.2 The per-reaction balance check

Each reaction's relative mass error is

$$
\epsilon_r = \frac{\big|\sum_j \texttt{product\_generated}_{j,r} - \sum_j \texttt{reactant\_consumed}_{j,r}\big|}
{\sum_j \texttt{reactant\_consumed}_{j,r}},
$$

with an absolute fallback when the denominator is $\le\texttt{eps()}$. Reaction
$r$ is balanced iff $\epsilon_r<\texttt{tol}$, default $\texttt{tol}=10^{-10}$.
Because the yields $\nu_{r,j}$ are validated at build time to sum to exactly $1.0$
(§6.8), this error should be at roundoff: $\sum_j\nu_{r,j}m_{\mathrm{flow},r}=m_{\mathrm{flow},r}=\sum_j\texttt{reactant\_consumed}$.
A nonzero $\epsilon_r$ would indicate a stoichiometry table that was not
normalized — exactly the bug class this independent ledger is meant to surface.

> **Implementation note.** `check_reaction_mass_balance` /
> `conservation_checks.jl:690–711`; `reaction_mass_flow_report`
> (lines 726–773) renders per-reaction and per-component breakdowns.

### 16.8.3 Scope and current limitation

The per-reaction ledger validates **stoichiometric** mass closure independently
of the global ledger. It does *not* yet feed reaction-specific heat release into
the energy ledger, so the energetic stoichiometry (whether the heat sources are
consistent with the reaction extents) is not independently checked — see §16.9.

---

## 16.9 What the discrete ledger does *not* yet accumulate

The discrete-form residual of §16.5 closes exactly only under *constant heat
capacity and constant cell volume*. Pyrolysis violates both: the mixture
composition $\xi$ evolves (changing $\rho c_p^\mathrm{eff}$ even at fixed $T$),
and in ALE mode the cell volumes $V_i$ change. Differentiating
$E_\mathrm{matrix}=\sum_i \rho c_{p,i}\,T_i\,V_i$ in time gives three terms:

$$
\frac{dE_\mathrm{matrix}}{dt}
= \sum_i \rho c_{p,i}\,V_i\,\frac{dT_i}{dt}
\;+\; \sum_i T_i\,V_i\,\frac{d(\rho c_{p,i})}{dt}
\;+\; \sum_i \rho c_{p,i}\,T_i\,\frac{dV_i}{dt} .
$$

Only the **first** term is what the source accumulators reconstruct (it is the
RHS times $\rho c_p V$). The matrix equation evolves $T$, *not* $\rho c_p\,T$, so
the discrete ledger is **missing**:

1. **Composition-change proxy term** $\displaystyle\int T\,\frac{\partial(\rho c_p^\mathrm{eff})}{\partial t}\,dV\,dt$
   — the change of the tracked $\rho c_p\,T$ proxy density at constant temperature as
   reactions convert high-$c_p$ virgin solid into low-$c_p$ char and release gas. (This
   term belongs to the ledger's proxy $E_\mathrm{matrix}$, not to a physical enthalpy
   density — the conserved quantity of the continuous formulation is
   $\sum_j\xi_j h_j(T)$, §3.3.2.) It is generally nonzero whenever reactions are
   active.
2. **ALE volume-change term** $\displaystyle\int \rho c_p^\mathrm{eff}\,T\,\frac{\partial V}{\partial t}\,dt$
   — the stored energy carried by a shrinking/expanding cell as its volume
   changes under mesh motion (§11).

Consequently the discrete-form residual as currently coded **cannot close to
integrator tolerance during active pyrolysis with time-varying composition or a
moving mesh**. Closing it would require per-substep, per-cell accumulation of the
*incremental* stored energy $\rho c_{p,i}(t)\,\Delta T_i\,V_i(t)$ (which absorbs
the composition and volume effects into the measured $\Delta(\rho c_p T)$), i.e.
wiring `record_substep!` into `compute_rhs!` so the tracker sees every cell's
state at every substep. That plumbing **does not exist in the package today**;
the `energy_ledger_audit.jl` example documents this explicitly and falls back to
the cleanly checkable quantities: the per-time surface energy-balance closure,
the gas-storage gap fraction, and the offline mass balance
(`examples/inverse_analysis/energy_ledger_audit.jl:1–37`).

This is the most important caveat for any user reading a runtime discrete
residual: a small *but nonzero* discrete residual during active reaction is
**expected**, not a bug, and reflects the unaccumulated composition/volume terms
— not an error in the flux/source assembly.

---

## 16.10 Callback infrastructure and runtime logging

### 16.10.1 DiagnosticsState and DiagnosticsRow

The running state of all three trackers, plus the emitted log, is held in a
single mutable container:

```julia
mutable struct DiagnosticsState{NC,NR}
    mass::MassConservationTracker{NC}
    energy::EnergyConservationTracker
    reactions::ReactionMassFlowTracker{NC,NR}
    stride::Int
    step_counter::Int
    last_t::Float64
    log::Vector{DiagnosticsRow}
end
```

It is constructed with `DiagnosticsState(mesh, material; stride=10)`, which seeds
each tracker from the current mesh state — establishing the $t=0$ baseline
($E_\mathrm{matrix}(0)$, $E_\mathrm{total}(0)$, $m_{\mathrm{initial},j}$)
(`callback.jl:47–77`). The default `stride` of 10 means one log row is emitted
per ten accepted substeps, trading output volume for temporal resolution.

Each log row is a concretely typed `NamedTuple`,

$$
\texttt{DiagnosticsRow} = \big(t,\ E_\mathrm{matrix},\ E_\mathrm{total},\
r_\mathrm{discrete},\ r_\mathrm{discrete}^\mathrm{rel},\
r_\mathrm{physical},\ r_\mathrm{physical}^\mathrm{rel},\
\Delta E_\mathrm{gas},\ \epsilon_{m}^{\max}\big),
$$

with all-`Float64` fields (`callback.jl:21–34`). Using a concrete `NamedTuple`
(rather than `Dict{Symbol,Any}`) keeps the log type-stable and means adding a new
diagnostic is a field addition, not a stringly-typed key.

### 16.10.2 record_substep! — the full cumulative path

`record_substep!` (`callback.jl:90–145`) is the entry point for the *full*
cumulative ledger. Per call it:

1. computes the substep size $\Delta t = t - \texttt{state.last\_t}$, updates
   `last_t`, and increments `step_counter`;
2. **guards** against $\Delta t\le0$ (a callback fired at a duplicate time, e.g.
   the initial point) by returning early — this prevents zero/negative-width
   substeps from polluting the accumulators;
3. calls `update_mass_tracker!` and `update_energy_tracker!` with the substep's
   boundary fluxes, source integrals, and area;
4. **only on stride boundaries** (`step_counter % stride == 0`) calls
   `check_energy_conservation` and `check_mass_conservation` against the current
   mesh and pushes a fully populated `DiagnosticsRow` (with real residuals).

This path requires the caller to supply, every substep, the conducted/applied
boundary fluxes, the three volume-integrated source totals, the transpiration
flux, the boundary mass-flux matrix, and the area — i.e. it must be threaded
through the residual assembly. As noted in §16.9, that threading is not yet
wired, so this path is exercised by tests rather than by `solve()`.

### 16.10.3 snapshot_diagnostics_row — the production path

What `solve()` *actually* uses is the **snapshot** path. `snapshot_diagnostics_row(t, mesh, material)`
(`callback.jl:157–171`) computes only the *instantaneous* quantities —
$E_\mathrm{matrix}(t)$, $E_\mathrm{total}(t)$, and $\Delta E_\mathrm{gas}(t)$ —
and fills **all cumulative fields with `NaN`**:

$$
r_\mathrm{discrete}=r_\mathrm{discrete}^\mathrm{rel}
=r_\mathrm{physical}=r_\mathrm{physical}^\mathrm{rel}
=\epsilon_m^{\max}=\texttt{NaN}.
$$

The `NaN` sentinel is a deliberate signal that those ledger residuals are *not
available in snapshot mode* (they require per-substep accumulation). The snapshot
path is cheap — three reductions over the mesh — and gives the user the two
quantities that *are* well-defined point-by-point: the energy content and the
gas-storage gap.

### 16.10.4 DiagnosticsCallback and solver integration

`DiagnosticsCallback(record!)` (`callback.jl:184–188`) wraps a user closure in a
`DiffEqCallbacks.DiscreteCallback` whose condition is `true` (fire on every
accepted step) and whose `affect!` runs `record!(integrator)`, with
`save_positions = (false, false)` so the callback does not force extra solution
saves.

In `solve()`, enabling diagnostics (`diagnostics=true`, default `false`;
`diagnostics_stride` default `10`, validated $\ge1$) installs a snapshot-mode
closure (`src/Solver/solve.jl:863–896`). The closure maintains its own
stride counter; on stride boundaries it (i) synchronizes the mesh from the
integrator's state with `update_mesh_from_state!(mesh, integrator.u, layout)` so
the trackers read the current cell states, then (ii) pushes a
`snapshot_diagnostics_row` to a local `Vector{DiagnosticsRow}`. After the solve,
the rows are appended to `solution.extras.diagnostics_log`
(`src/Solver/solve.jl:937–939`; the `extras` schema reserves
`diagnostics_log` at `src/Problem/pyrolysis_problem.jl:275`). Reading
`sol.extras.diagnostics_log` therefore yields the energy and gas-gap time series,
with `NaN` in the cumulative-residual columns under the current production wiring.

---

## 16.11 Comparison to ThermaKin2Ds, Gpyro, and FDS

The conservation-diagnostics design follows the **FDS/Gpyro family** in its
modeling convention but adds an explicit, runtime, dual-ledger audit that is not
a standard feature of any of the three reference codes.

**Heat-capacity / gas-storage convention.** Pyrolysis.jl excludes gas from
$\rho c_p^\mathrm{eff}$ and carries gas enthalpy by advection — the **FDS and
Gpyro convention**. This differs from **ThermaKin2Ds**, whose energy equation
(ThermaKin Eq. 17, mapping their summed heat capacity to our
$\sum_j \xi_j c_{p,j}$) *keeps* gas storage in the volumetric heat capacity (all
components summed). Pyrolysis.jl's physical-form ledger is, in effect, a runtime
measurement of the energy difference between the FDS/Gpyro convention it uses and
the ThermaKin convention it does not — the gas-storage gap $\Delta E_\mathrm{gas}$
is precisely the term ThermaKin retains and Pyrolysis.jl drops.

**Energy accounting.** **FDS** performs an internal solid-phase energy budget and
exposes global energy-conservation accounting (the FDS Technical Reference
Guide's energy-budget appendix, Appendix P analog), reporting cumulative
boundary, radiative, and chemical contributions against stored enthalpy.
Pyrolysis.jl's *discrete-form ledger* is the direct analog at the cell level:
both reconstruct stored energy from accumulated boundary and source fluxes and
report the closure residual. Where FDS folds everything into one budget,
Pyrolysis.jl deliberately splits the budget into the discrete (bookkeeping) and
physical (model-cost) ledgers so the two failure modes — assembly bugs vs.
approximation breakdown — are distinguishable.

**Mass accounting.** **Gpyro** advances control-volume conservation of
$\bar\rho\,\Delta z$ and $\bar\rho Y_i\,\Delta z$ with a fully-implicit
backward-Euler/Patankar scheme; conservation is enforced by the discretization
itself, with closure verified through the solid-density bookkeeping. Pyrolysis.jl
instead conserves species mass by construction (validated stoichiometry, §6.8)
and *audits* it with the independent global and per-reaction ledgers. The
per-reaction ledger has no direct counterpart in the reference codes; it is a
build-time-plus-runtime guard specific to this code's single-reactant,
yield-normalized stoichiometry model.

**Reference-point deformation vs. ALE.** ThermaKin handles condensed-phase
deformation through a single-reference-point density-integral source term, and
Gpyro through a variable-$\Delta z$ control volume; Pyrolysis.jl uses ALE mesh
motion (§10, §11). This is why the Pyrolysis.jl ledgers must carry a
time-dependent area $A(t)$ (trapezoidal average over a substep) and why the
unaccumulated ALE volume-change term (§16.9) is a closure caveat unique to the
ALE formulation that the reference codes, with their different deformation
treatments, do not incur in the same form.

---

## 16.12 Worked summary of the ledger equations

For quick reference, the complete set of accumulated and checked quantities:

| Ledger | Stored quantity | Expected | Residual | Tolerance |
|---|---|---|---|---|
| Mass (per comp. $j$) | $m_j=\sum_{i\in\text{act}}\xi_{j,i}V_i$ | $m_{\mathrm{init},j}-m_{\mathrm{efflux},j}$ | $m_j-m_{\mathrm{init},j}+m_{\mathrm{efflux},j}$ | $10^{-10}$ (rel.; abs. if $m_\mathrm{init}\approx0$) |
| Energy, discrete | $E_\mathrm{matrix}=\sum \rho c_p^\mathrm{eff}T_iV_i$ | $E_\mathrm{m}(0)+\!\int\!(q_\mathrm{cond}A+Q_\mathrm{rxn}+S_\mathrm{conv}+S_\mathrm{gen})$ | $E_\mathrm{matrix}-E_\mathrm{matrix}^\mathrm{exp}$ | $10^{-6}$ (rel.) |
| Energy, physical | $E_\mathrm{total}=\sum \rho c_p^\mathrm{total}T_iV_i$ | $E_\mathrm{t}(0)+\!\int\!(q_\mathrm{BC}A+Q_\mathrm{rxn})-\!\int\! q_\mathrm{transp}A$ | $E_\mathrm{total}-E_\mathrm{total}^\mathrm{exp}$ | $<1\%$ (canonical) |
| Reaction (per rxn $r$) | $\sum_j$ generated, $\sum_j$ consumed | generated $=$ consumed | $|\text{gen}-\text{cons}|/\text{cons}$ | $10^{-10}$ (rel.) |
| Gas gap | $\Delta E_\mathrm{gas}=E_\mathrm{total}-E_\mathrm{matrix}$ | — (diagnostic) | — | $\ll1\%$ of $E_\mathrm{matrix}$ |

Sign-convention reminders carried through every row: $J>0$ = inflow (so efflux
$=-J A\Delta t$); divergence positive = outflow; $Q_\mathrm{rxn}=-\sum_r h_r r_r$
with $h>0$ endothermic; $S_\mathrm{conv}$ typically negative; transpiration enters
$E^\mathrm{exp}_\mathrm{total}$ with a minus sign (energy leaving).

---

## 16.13 Limitations and open issues

The following are the known gaps in the conservation diagnostics, drawn from the
source documentation, the audit example, and the KB open questions. They are
caveats on *interpretation*, not correctness bugs in the kernels.

1. **Composition-varying discrete closure (primary limitation).** As derived in
   §16.9, the discrete-form residual omits the composition-change proxy term
   $\int T\,\partial(\rho c_p^\mathrm{eff})/\partial t\,dV\,dt$ and the ALE
   volume-change term $\int \rho c_p^\mathrm{eff}\,T\,\partial V/\partial t\,dt$.
   A residual that closes to integrator tolerance during active pyrolysis
   requires per-substep, per-cell accumulation of $\rho c_{p,i}\,\Delta T_i\,V_i$,
   i.e. wiring `record_substep!` into `compute_rhs!`. This is **not implemented**;
   in production, `solve()` runs the snapshot path with `NaN` residual columns.
   The audit example
   (`examples/inverse_analysis/energy_ledger_audit.jl`) documents the gap and
   substitutes per-time surface-balance, gas-gap, and offline-mass checks.

2. **Snapshot incompleteness.** Because `solve()` uses
   `snapshot_diagnostics_row`, the production `diagnostics_log` carries only
   $E_\mathrm{matrix}$, $E_\mathrm{total}$, and $\Delta E_\mathrm{gas}$; the
   discrete/physical residuals and the mass-error column are `NaN`. The full
   cumulative ledger is available only by manually threading `record_substep!`.

3. **Per-reaction energetic closure.** The per-reaction ledger validates
   stoichiometric mass closure but does not yet integrate reaction-specific heat
   release into the energy ledger, so energetic stoichiometry is not
   independently audited (KB open question; §16.8.3).

4. **Boundary-flux decomposition.** The energy ledger treats $q_\mathrm{BC,applied}$
   as a single scalar per boundary. Decomposing it into separate radiation,
   convection, and conduction channels would let each mechanism be audited
   independently — the solver already exposes the *components* in `extras`
   (`surface_heat_flux_absorbed`, `_reradiation`, `_convection`, `_conduction`,
   `_transpiration`), but the tracker consumes only the aggregate.

5. **Species-specific transpiration.** Transpiration is summed over all gas
   species in the energy ledger. Tracking per-species transpiration enthalpy
   would expose if one component dominates the escape and would let the
   multi-component diffusion model be validated energetically (KB open question).

6. **Asymmetric ambient temperatures with a permeable back face.** The
   generation-sink/transpiration telescoping (§16.6.3) uses a single scalar
   $T_\mathrm{ref}$ resolved from one boundary's ambient. If the two boundaries
   have different ambients *and* the back face is permeable, the ledger fails to
   close at the bottom by $\sum_g \dot m_{g,\mathrm{bottom}}c_{p,g}(T_{\infty,\mathrm{top}}-T_{\infty,\mathrm{bottom}})$.
   This is harmless for the standard impermeable back face (both terms vanish) but
   would require a per-cell datum to support in general (see §7.5 caveat).

7. **Pressure scaling of the gas gap.** The gas-storage gap (and hence the
   physical residual) scales with pressure; no automatic switch to the full gas
   storage at elevated pressure is implemented. The diagnostic makes the cost
   *visible* (the user is warned when $r_\mathrm{physical}>1\%$) but does not
   *act* on it.

---

*Cross-references.* The governing equations these ledgers audit are stated in §3
(Governing Equations); the energy-source terms $S_\mathrm{conv}$, $S_\mathrm{gen}$,
the matrix/total heat capacities, and the $S_\mathrm{conv}$/transpiration mutual
exclusivity in §7 (Heat Conduction, Convection & Energy Assembly); the
heat-of-reaction sign and stoichiometric closure in §6 (Reaction Kinetics); the
surface energy balance and boundary fluxes in §12 (Boundary & Initial
Conditions); the time-dependent area and volume from §10 (Volume Change) and §11
(ALE Framework); the active-cell mask from §14 (Adaptive Mesh Refinement); and
the callback/integration machinery from §15 (Temporal Integration). The
diagnostics are themselves a verification tool used in §18 (Verification &
Validation).
