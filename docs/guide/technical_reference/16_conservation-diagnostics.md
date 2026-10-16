# 16. Conservation Diagnostics

This chapter derives the runtime conservation diagnostics that Pyrolysis.jl uses
to audit its own numerics. We develop a per-component **mass ledger** and a
**sensible-energy ledger** that closes to integrator tolerance by construction
for constant composition and volume. For every accumulator we give the
continuous balance law, its exact discrete counterpart as the code assembles it,
the precise sign conventions, and the conditions under which the residual is
expected to vanish. The chapter closes by stating exactly which balance terms
are *not* yet accumulated, why the energy residual cannot close for time-varying
composition without them, and how the production diagnostics callback operates
in snapshot mode.

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
transport, a moving (ALE) mesh, and surface-cell depletion, all advanced by an
implicit time integrator. With this many interacting mechanisms, a silent error
in flux or source assembly can masquerade as physics. The diagnostics subsystem
exists to make conservation **observable at runtime** rather than only in offline
post-processing, so that a broken term announces itself as a non-closing ledger.

The question the energy ledger answers is one of bookkeeping: do the
accumulated boundary, reaction, and gas-coupling terms add up to the change in
stored energy that the RHS actually evolves? It is designed to close to
integrator tolerance *by construction* — it is assembled from exactly the terms
the RHS computes. A failure is a bug in flux/source assembly or in the callback
plumbing, not a physics error.

The key design subtlety, stated plainly in the source documentation of
`EnergyConservationTracker`, is that this ledger audits *bookkeeping, not the
physical model*: because it is built from the terms the RHS actually assembles,
**it will close even if the RHS itself double-counts a mechanism.** Catching that
class of error is the job of the independent closure tests (the sealed-column and
gas-enthalpy verification cases of §18) and of the formulation itself (the surface
balance carries no gas-enthalpy term alongside $S_\mathrm{conv}$, §7.4.3).

> **Implementation note.** The two trackers and the callback live in module
> `Diagnostics`. They are pure functions of the *current mesh state*
> (`mesh.cell_states`, `mesh.cell_volumes`, `mesh.active_cell_mask`) plus a small
> set of mutable accumulators; they never read the raw state vector directly. The
> solver synchronizes the mesh from the integrator's `u` before each diagnostic
> snapshot via `update_mesh_from_state!` (`src/Solver/solve.jl`).

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
build-time stoichiometric closure (§6.5.3) to guarantee $\sum_j S_j=0$. Under ALE
every component rides the mesh (§11.1.2), so $J_j$ — the flux relative to the
skeleton — is the only transport channel on a moving mesh too. With the two
boundaries (top $z=L$, bottom $z=0$) and the sign convention that $J_j>0$ is
$+z$ transport, the *outward* normal flux is $+J_j$ at the top and $-J_j$ at the
bottom. The diagnostics layer adopts the **into-domain** convention for its
boundary inputs (positive $J$ = inflow), under which the net efflux is
$-(J_{\mathrm{top},j}+J_{\mathrm{bot},j})A$ — see §16.3.

### 16.2.2 Energy

Multiplying the temperature equation (§7.3) by $\rho c_p^\mathrm{eff}\,V$ and
integrating, the sensible energy the equation evolves obeys

$$
\frac{d}{dt}\int_\Omega \rho c_p^\mathrm{eff}\,T\,dV
\;\approx\; \int_{\partial\Omega} q_\mathrm{cond}\,dA
+ \int_\Omega Q_\mathrm{rxn}\,dV
+ \int_\Omega S_\mathrm{conv}\,dV ,
$$

with $\rho c_p^\mathrm{eff}=\sum_j \xi_j c_{p,j}(T)$ the mixture heat capacity over
**every** component, gas included (§7.3.2), and $q_\mathrm{cond}$ the *conducted*
boundary flux. The surface balance has no gas-enthalpy term, so the conducted
flux equals the applied boundary flux: the gas crosses the surface at $T_s$ on
both sides (§1.4.3). The volumetric radiation source $\int_\Omega Q_\mathrm{rad}\,dV$ of
the temperature equation is left out (it is zero except under `BEER_LAMBERT`;
§16.11). This is the law the energy ledger
discretizes. The $\approx$ flags two terms missing from the storage
time-derivative — the composition-change term and the ALE volume-change term —
which are the subject of §16.7.

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
(`mesh.active_cell_mask[i]`) excludes parked cells (see §14, Adaptive Mesh
Refinement), so the ledger always reflects only the live solution.

> **Implementation note.** `compute_domain_mass`,
> `src/Diagnostics/conservation_checks.jl`. The loop iterates
> `eachindex(mesh.cells)` and `continue`s on inactive cells; it accumulates into
> a `MVector{NC,Float64}` and returns an `NTuple{NC,Float64}`. The mask test must
> be performed cell-by-cell — never `eachindex` over a compacted array — which is
> the AMR invariant noted in §14.

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

> **Implementation note.** `update_mass_tracker!`,
> `conservation_checks.jl`. `J_boundary` is an `NC×2` matrix with column 1
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
equals `mesh.cross_section_area_initial`.

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

**Scope of the $10^{-10}$ tolerance.** The roundoff-level closure holds wherever
the discrete ledger telescopes exactly: every component on a fixed (Eulerian)
mesh, and every component under ALE, where the only mesh-motion term is the
dilation $-\xi_j\theta_i$ (which conserves each cell's mass $\xi_jV$) and the
transport flux $J_j$ at the boundary faces is the whole of the efflux (§11.2.3).
A sealed collapsing column conserves each component to integrator tolerance
(V-CPL-03 P4).

> **Implementation note.** `check_mass_conservation` and
> `mass_conservation_report`, `conservation_checks.jl`; the latter renders the
> per-component table and an overall PASSED/FAILED line. The state is held in
> `MassConservationTracker{NC}`: `initial_mass`, `cumulative_efflux`,
> `last_time`.

### 16.3.4 Sources of (small) mass-ledger discrepancy

Even when the implementation is correct the ledger carries avoidable error
sources. (i) **Offline vs. online integration.** If the efflux is reconstructed
offline by trapezoidal integration of saved surface fluxes (rather than
accumulated online per substep), the integration error of the post-processing
quadrature appears in the residual. (ii) **Substep granularity.** The online
accumulator uses the substep-mean area; on a rapidly shrinking ALE mesh the
trapezoidal-area approximation is second-order in $\Delta t$. (iii) **Active-mask
transitions.** When a surface cell is parked (§14.3), its mass is merged into
its neighbour exactly; the ledger picks up the post-merge masses, so a correctly
conserving merge leaves no residual but a buggy one would surface here.

---

## 16.4 Energy storage

The stored quantity is the sensible energy the temperature equation evolves,

$$
E(t) = \sum_{i\in\mathrm{active}} \rho c_p^\mathrm{eff}(\xi_i, T_i)\, T_i\, V_i,
\qquad
\rho c_p^\mathrm{eff}(\xi,T)=\sum_{j=1}^{N_C}\xi_j\,c_{p,j}(T),
$$

with the mixture heat capacity over every component (`@generated` summation
kernel in `effective_properties.jl`, no phase branch). The sum respects the
active-cell mask. The product $\rho c_p\cdot T$ uses the **absolute** temperature
$T$ (not $T-T_\mathrm{ref}$): the diagnostics measure a sensible-energy *content*
relative to $0\,\mathrm{K}$, which is acceptable because the ledger only ever uses
*differences* $E(t)-E(0)$ in which the datum cancels (provided $\rho c_p$ is held
fixed; see the composition-change caveat in §16.7).

> **Implementation note.** `compute_sensible_energy`,
> `conservation_checks.jl`, calling `Properties.effective_heat_capacity`.

---

## 16.5 The energy ledger

### 16.5.1 Continuous-to-discrete construction

The ledger is assembled from *exactly the terms the RHS computes*. Integrating
the energy balance of §16.2.2 in time and discretizing each term by the same
quadrature the RHS uses gives

$$
E(t) = E(0)
+ \int_0^t \Big[\, q_\mathrm{cond,boundary}\,A
+ Q_\mathrm{rxn}^\mathrm{tot}
+ S_\mathrm{conv}^\mathrm{tot} \,\Big]\,d\tau,
$$

where each integrand is the *volume-integral of the corresponding RHS field*:

- $q_\mathrm{cond,boundary}$ — net **conducted** boundary heat flux [W·m⁻²],
  positive into the slab. This is the flux the Newton boundary solver drives
  into the first/last cell (§12.2); the surface balance has no gas-enthalpy
  term, so it equals the applied BC flux.
- $Q_\mathrm{rxn}^\mathrm{tot}=\sum_i Q_{\mathrm{rxn},i}V_i = \int_\Omega Q_\mathrm{rxn}\,dV$ [W],
  the total reaction heat ($Q_\mathrm{rxn}=-\sum_r h_r r_r$; §6.6).
- $S_\mathrm{conv}^\mathrm{tot}=\int_\Omega S_\mathrm{conv}\,dV$ [W], the total gas
  advective enthalpy sink; **typically negative** (heat carried out by gas; §7.4).

### 16.5.2 Per-substep accumulation

The energy tracker accumulates the time-integrals as Euler products over each
substep $\Delta t$:

$$
\begin{aligned}
\texttt{cum\_boundary\_conduction} &\mathrel{+}= q_\mathrm{cond,boundary}\,A\,\Delta t,\\
\texttt{cum\_reaction\_heat}       &\mathrel{+}= Q_\mathrm{rxn}^\mathrm{tot}\,\Delta t,\\
\texttt{cum\_advection\_sink}      &\mathrel{+}= S_\mathrm{conv}^\mathrm{tot}\,\Delta t .
\end{aligned}
$$

The boundary term carries the area factor $A$ and has flux units [W·m⁻²], so
$q\cdot A\cdot\Delta t$ has units of joules; the source terms are already
volume-integrals [W], so $\cdot\,\Delta t$ gives joules.

> **Implementation note.** `update_energy_tracker!`, `conservation_checks.jl`.
> The accumulators live in `EnergyConservationTracker`: `initial_energy`,
> `cumulative_boundary_conduction`, `cumulative_reaction_heat`,
> `cumulative_advection_sink`, `last_time`.

### 16.5.3 The residual

At a check point the *expected* energy is the initial energy plus the
accumulated sources,

$$
E^\mathrm{exp}(t) = E(0)
+ \texttt{cum\_boundary\_conduction}
+ \texttt{cum\_reaction\_heat}
+ \texttt{cum\_advection\_sink},
$$

and the **residual** is the gap between the measured and expected stored energy:

$$
\boxed{\,r(t) = E(t) - E^\mathrm{exp}(t)\,}
$$

with relative magnitude $r/|E(0)|$ (absolute if $|E(0)|\le\texttt{eps()}$) and
the verdict $r^\mathrm{rel}<\texttt{tol}$, default $\texttt{tol}=10^{-6}$.

> **Implementation note.** `check_energy_conservation`, `conservation_checks.jl`,
> returning `residual`, `relative`, `conserved`, `current_energy`.

### 16.5.4 Why it closes — and what it does *not* prove

This residual closes to integrator tolerance *by construction* whenever the
heat capacity and volume are constant over the trajectory and no radiation is
absorbed in depth, because both
$E(t)$ and its expectation are assembled from the same RHS fields under exactly
consistent quadrature; the only error is the integrator's local truncation error
in advancing $T$, which the implicit solver controls to `reltol`. A non-closing
residual therefore points to a **bookkeeping bug**: a flux computed but not
accumulated, a callback firing at the wrong time, a sign error in a source term,
or a missing factor of $A$ or $V$.

The crucial limitation, repeated from §16.1: the ledger **cannot detect a
double-count**. If the RHS were (incorrectly) to add a surface gas-enthalpy term
on top of $S_\mathrm{conv}$, the accumulator would faithfully record both, the
expectation would include both, and the residual would still close. Catching that
class of error is the job of the independent verification cases (§18). The ledger is also subject to the
composition-change and ALE volume-change omissions of §16.7, which prevent
*exact* closure whenever the mixture heat capacity or cell volume vary in time.

---

## 16.6 Boundary-flux semantics

The energy ledger accumulates the **conducted flux** $q_\mathrm{cond,boundary}$:
the heat actually conducted into the first/last interior cell, i.e.
$k_\mathrm{eff}(T_s-T_\mathrm{cell})/\Delta z_{1/2}$ evaluated at the converged
surface temperature $T_s$ from the Newton boundary solve (§12.2). This is what the
energy equation receives. The surface energy balance the Newton solver enforces
is (§12.2)

$$
q_\mathrm{BC,applied}(T_s) = q_\mathrm{cond}(T_s),
$$

with no gas-enthalpy term (§1.4.3, §7.4.3), so the conducted flux *is* the
applied flux: the externally
driven flux that the boundary condition imposes (for a `RadiativeFluxBC` the
absorbed incident radiation less the re-radiation and convective exchange; for a
pure `HeatFluxBC` the prescribed $q(t)$). The solver exposes these as
`extras.surface_heat_flux_net` ($q_\mathrm{BC,applied}$ — the applied face flux in
every radiation regime; under `BEER_LAMBERT` it excludes the incident radiation,
which enters volumetrically), `extras.surface_heat_flux_conduction`
($q_\mathrm{cond}$), plus `extras.volumetric_absorption` (the Beer–Lambert in-depth
deposit $\alpha_{\text{surf}}q_{\text{ext}}(1-e^{-\tau_{\text{tot}}})$; zero for other
models) and `extras.back_heat_flux` for the bottom boundary
(`src/Problem/pyrolysis_problem.jl`).

---

## 16.7 What the energy ledger does *not* yet accumulate

The residual of §16.5 closes exactly only under *constant heat capacity and
constant cell volume*. Pyrolysis violates both: the mixture composition $\xi$
evolves (changing $\rho c_p^\mathrm{eff}$ even at fixed $T$), and in ALE mode the
cell volumes $V_i$ change. Differentiating $E=\sum_i \rho c_{p,i}\,T_i\,V_i$ in
time gives three terms:

$$
\frac{dE}{dt}
= \sum_i \rho c_{p,i}\,V_i\,\frac{dT_i}{dt}
\;+\; \sum_i T_i\,V_i\,\frac{d(\rho c_{p,i})}{dt}
\;+\; \sum_i \rho c_{p,i}\,T_i\,\frac{dV_i}{dt} .
$$

Only the **first** term is what the source accumulators reconstruct (it is the
RHS times $\rho c_p V$). The temperature equation evolves $T$, *not* $\rho c_p\,T$,
so the ledger is **missing**:

1. **Composition-change proxy term** $\displaystyle\int T\,\frac{\partial(\rho c_p^\mathrm{eff})}{\partial t}\,dV\,dt$
   — the change of the tracked $\rho c_p\,T$ proxy density at constant temperature as
   reactions convert high-$c_p$ virgin solid into low-$c_p$ char and release gas,
   and as the gas concentration changes by transport and dilation. (This term
   belongs to the ledger's proxy $E$, not to a physical enthalpy density — the
   conserved quantity of the continuous formulation is $\sum_j\xi_j h_j(T)$,
   §3.3.2.) It is generally nonzero whenever reactions are active.
2. **ALE volume-change term** $\displaystyle\int \rho c_p^\mathrm{eff}\,T\,\frac{\partial V}{\partial t}\,dt$
   — the stored energy carried by a shrinking/expanding cell as its volume
   changes under mesh motion (§11). For the dilation part this cancels the
   corresponding piece of term 1 exactly (dilation conserves $\xi_jV$).

Consequently the residual as currently coded **cannot close to integrator
tolerance during active pyrolysis with time-varying composition or a moving
mesh**. Closing it would require per-substep, per-cell accumulation of the
*incremental* stored energy $\rho c_{p,i}(t)\,\Delta T_i\,V_i(t)$ (which absorbs
the composition and volume effects into the measured $\Delta(\rho c_p T)$), i.e.
wiring `record_substep!` into `compute_rhs!` so the tracker sees every cell's
state at every substep. That plumbing **does not exist in the package**; the
verification case V-CPL-06 carries its own per-step ledger with the hand-off and
dilation terms made explicit and closes the ALE energy budget to $5\times10^{-5}$
of the boundary input (§18.4).

This is the most important caveat for any user reading a runtime residual: a
small *but nonzero* residual during active reaction is **expected**, not a bug,
and reflects the unaccumulated composition/volume terms — not an error in the
flux/source assembly.

---

## 16.8 Callback infrastructure and runtime logging

### 16.8.1 DiagnosticsState and DiagnosticsRow

The running state of both trackers, plus the emitted log, is held in a
single mutable container:

```julia
mutable struct DiagnosticsState{NC}
    mass::MassConservationTracker{NC}
    energy::EnergyConservationTracker
    stride::Int
    step_counter::Int
    last_t::Float64
    log::Vector{DiagnosticsRow}
end
```

It is constructed with `DiagnosticsState(mesh, material; stride=10)`, which seeds
each tracker from the current mesh state — establishing the $t=0$ baseline
($E(0)$, $m_{\mathrm{initial},j}$). The default `stride` of 10 means one log row
is emitted per ten accepted substeps, trading output volume for temporal
resolution.

Each log row is a concretely typed `NamedTuple`,

$$
\texttt{DiagnosticsRow} = \big(t,\ E_\mathrm{sensible},\ r_\mathrm{energy},\ r_\mathrm{energy}^\mathrm{rel},\ \epsilon_{m}^{\max}\big),
$$

with all-`Float64` fields (`callback.jl`). Using a concrete `NamedTuple` (rather
than `Dict{Symbol,Any}`) keeps the log type-stable and means adding a new
diagnostic is a field addition, not a stringly-typed key.

### 16.8.2 record_substep! — the full cumulative path

`record_substep!` (`callback.jl`) is the entry point for the *full* cumulative
ledger. Per call it:

1. computes the substep size $\Delta t = t - \texttt{state.last\_t}$ and updates
   `last_t`;
2. **guards** against $\Delta t\le0$ (a callback fired at a duplicate time, e.g.
   the initial point) by returning early, before `step_counter` is incremented — a
   duplicate time neither pollutes the accumulators nor advances the stride;
3. increments `step_counter` and calls `update_mass_tracker!` and `update_energy_tracker!` with the substep's
   conducted boundary flux, the two volume-integrated source totals
   ($Q_\mathrm{rxn}$, $S_\mathrm{conv}$), the boundary mass-flux matrix, and the
   area;
4. **only on stride boundaries** (`step_counter % stride == 0`) calls
   `check_energy_conservation` and `check_mass_conservation` against the current
   mesh and pushes a fully populated `DiagnosticsRow` (with real residuals).

This path requires the caller to supply, every substep, the conducted boundary
flux, the two volume-integrated source totals, the boundary mass-flux matrix,
and the area — i.e. it must be threaded through the residual assembly. As noted
in §16.7, that threading is not yet wired, so this path is exercised by tests
rather than by `solve()`.

### 16.8.3 snapshot_diagnostics_row — the production path

What `solve()` *actually* uses is the **snapshot** path.
`snapshot_diagnostics_row(t, mesh, material)` (`callback.jl`) computes only the
*instantaneous* sensible energy $E(t)$ and fills the cumulative fields with
`NaN`:

$$
r_\mathrm{energy}=r_\mathrm{energy}^\mathrm{rel}=\epsilon_m^{\max}=\texttt{NaN}.
$$

The `NaN` sentinel is a deliberate signal that those ledger residuals are *not
available in snapshot mode* (they require per-substep accumulation). The snapshot
path is cheap — one reduction over the mesh — and gives the user the quantity
that *is* well-defined point-by-point: the energy content.

### 16.8.4 DiagnosticsCallback and solver integration

`DiagnosticsCallback(record!)` (src/Diagnostics/callback.jl) wraps a user
closure in a `SciMLBase.DiscreteCallback` whose condition is `true` (fire on
every accepted step) and whose `affect!` runs `record!(u, t)` on the
integrator's state and time, with `save_positions = (false, false)` so the
callback does not force extra solution saves. `record!` sits behind a
fixed-signature `FunctionWrapper`, so the callback's type does not depend on the
closure's (§15.5.3).

In `solve()`, enabling diagnostics (`diagnostics=true`, default `false`;
`diagnostics_stride` default `10`, validated $\ge1$) installs a snapshot-mode
closure (`solve`, src/Solver/solve.jl). The closure maintains its own
stride counter; on stride boundaries it (i) synchronizes the mesh from the
integrator's state with `update_mesh_from_state!(mesh, u, layout)` so
the trackers read the current cell states, then (ii) pushes a
`snapshot_diagnostics_row` to a local `Vector{DiagnosticsRow}`. After the solve,
the rows are appended to `solution.extras.diagnostics_log` (the `extras`
schema, `_extras` in src/Problem/pyrolysis_problem.jl, reserves
`diagnostics_log`). Reading `sol.extras.diagnostics_log` therefore yields the
energy time series, with `NaN` in the cumulative-residual columns under the
current production wiring.

---

## 16.9 Comparison to ThermaKin2Ds, Gpyro, and FDS

The conservation-diagnostics design adds an explicit, runtime ledger audit that
is not a standard feature of any of the three reference codes.

**Energy accounting.** **FDS** performs an internal solid-phase energy budget and
exposes global energy-conservation accounting (the FDS Technical Reference
Guide's energy-budget appendix, Appendix P analog), reporting cumulative
boundary, radiative, and chemical contributions against stored enthalpy.
Pyrolysis.jl's energy ledger is the direct analog at the cell level: both
reconstruct stored energy from accumulated boundary and source fluxes and report
the closure residual.

**Mass accounting.** **Gpyro** advances control-volume conservation of
$\bar\rho\,\Delta z$ and $\bar\rho Y_i\,\Delta z$ with a fully-implicit
backward-Euler/Patankar scheme; conservation is enforced by the discretization
itself, with closure verified through the solid-density bookkeeping. Pyrolysis.jl
instead conserves species mass by construction (stoichiometry validated when each
reaction is built, §6.5.3) and *audits* it with the global mass ledger.

**Reference-point deformation vs. ALE.** ThermaKin handles condensed-phase
deformation through a single-reference-point density-integral source term, and
Gpyro through a variable-$\Delta z$ control volume; Pyrolysis.jl uses ALE mesh
motion (§10, §11). This is why the Pyrolysis.jl ledgers must carry a
time-dependent area $A(t)$ (trapezoidal average over a substep) and why the
unaccumulated ALE volume-change term (§16.7) is a closure caveat unique to the
ALE formulation that the reference codes, with their different deformation
treatments, do not incur in the same form.

---

## 16.10 Worked summary of the ledger equations

For quick reference, the complete set of accumulated and checked quantities:

| Ledger | Stored quantity | Expected | Residual | Tolerance |
|---|---|---|---|---|
| Mass (per comp. $j$) | $m_j=\sum_{i\in\text{act}}\xi_{j,i}V_i$ | $m_{\mathrm{init},j}-m_{\mathrm{efflux},j}$ | $m_j-m_{\mathrm{init},j}+m_{\mathrm{efflux},j}$ | $10^{-10}$ (rel.; abs. if $m_\mathrm{init}\approx0$) |
| Energy | $E=\sum \rho c_p^\mathrm{eff}T_iV_i$ | $E(0)+\!\int\!(q_\mathrm{cond}A+Q_\mathrm{rxn}+S_\mathrm{conv})$ | $E-E^\mathrm{exp}$ | $10^{-6}$ (rel.) |

Sign-convention reminders carried through every row: $J>0$ = inflow (so efflux
$=-J A\Delta t$); divergence positive = outflow; $Q_\mathrm{rxn}=-\sum_r h_r r_r$
with $h>0$ endothermic; $S_\mathrm{conv}$ typically negative.

---

## 16.11 Limitations and open issues

The following are the known gaps in the conservation diagnostics. They are
caveats on *interpretation*, not correctness bugs in the kernels.

1. **Composition-varying closure (primary limitation).** As derived in §16.7,
   the residual omits the composition-change proxy term
   $\int T\,\partial(\rho c_p^\mathrm{eff})/\partial t\,dV\,dt$ and the ALE
   volume-change term $\int \rho c_p^\mathrm{eff}\,T\,\partial V/\partial t\,dt$.
   A residual that closes to integrator tolerance during active pyrolysis
   requires per-substep, per-cell accumulation of $\rho c_{p,i}\,\Delta T_i\,V_i$,
   i.e. wiring `record_substep!` into `compute_rhs!`. This is **not implemented**;
   in production, `solve()` runs the snapshot path with `NaN` residual columns.

2. **Snapshot incompleteness.** Because `solve()` uses
   `snapshot_diagnostics_row`, the production `diagnostics_log` carries only
   $E_\mathrm{sensible}$; the residual and the mass-error column are `NaN`. The
   full cumulative ledger is available only by manually threading
   `record_substep!`.

3. **Boundary-flux decomposition.** The energy ledger treats the boundary flux
   as a single scalar per boundary. Decomposing it into separate radiation,
   convection, and conduction channels would let each mechanism be audited
   independently — the solver already exposes the *components* in `extras`
   (`surface_heat_flux_absorbed`, `_reradiation`, `_convection`, `_conduction`),
   but the tracker consumes only the aggregate.

4. **No radiation term.** The energy ledger accumulates the conducted boundary
   flux, $Q_\mathrm{rxn}$, and $S_\mathrm{conv}$ only; a volumetric Beer–Lambert
   deposit $Q_\mathrm{rad}$ is not among its accumulators, so under `BEER_LAMBERT`
   the residual also carries the absorbed in-depth radiation.

---

*Cross-references.* The governing equations these ledgers audit are stated in §3
(Governing Equations); the energy-source term $S_\mathrm{conv}$, the mixture heat
capacity, and the absence of a surface gas-enthalpy term in §7 (Heat Conduction,
Convection & Energy Assembly); the
heat-of-reaction sign and stoichiometric closure in §6 (Reaction Kinetics); the
surface energy balance and boundary fluxes in §12 (Boundary & Initial
Conditions); the time-dependent area and volume from §10 (Volume Change) and §11
(ALE Framework); the active-cell mask from §14 (Adaptive Mesh Refinement and
Surface Depletion); and
the callback/integration machinery from §15 (Temporal Integration). The
diagnostics are themselves a verification tool used in §18 (Verification &
Validation).
