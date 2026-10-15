# 17. Sensitivity Analysis (Theory)

This chapter develops the mathematics of parametric sensitivity for the
Pyrolysis.jl solver: how the package computes the derivative of a simulation
output with respect to model parameters, why those derivatives are
well-defined, and how they support inverse parameter estimation, uncertainty
quantification, and gradient-based optimization. We treat both *local*
sensitivity (forward-mode algorithmic differentiation through the entire
`solve` chain, and reverse-mode adjoint differentiation of a scalar loss) and
*global*, variance-based sensitivity (Sobol indices computed through the
parameter-sweep harness). Throughout, the central enabling fact is that the
unified residual of §15 is written to be smooth and scalar-type-generic: the
state cache, the state-derived geometry, the boundary-temperature Newton solve,
and every property evaluation carry the scalar element type `Tu` end to end, so
a `ForwardDiff.Dual` (or a Zygote pullback) propagates cleanly from a parameter
all the way to a quantity of interest. The implementation lives in the optional
`PyrolysisSciMLSensitivityExt` extension behind two stub functions; the
mathematics below is what those wrappers actually differentiate.

Conventions used in this chapter (binding, from §2, Nomenclature):

- $\boldsymbol{\theta}\in\mathbb{R}^{p}$ is the vector of **parameters** with
  respect to which sensitivities are taken (e.g. $A_r$, $E_r$, $h_r$, a density
  $\rho_j$, a conductivity $k_j$, or a boundary-condition coefficient). $p$ is
  the parameter count.
- $\mathbf{u}(t)\in\mathbb{R}^{N_{\text{state}}}$ is the flat ODE state vector
  packed block-major as $[\,T \mid \xi_1\dots\xi_{N_C}\mid z\mid\chi\,]$ (§15.1).
- $\mathbf{y}=\mathcal{O}(\mathbf{u},\boldsymbol{\theta})\in\mathbb{R}^{m}$ is a
  vector **output** (a quantity of interest, "QoI") extracted from the
  solution; $\ell(\boldsymbol{\theta})\in\mathbb{R}$ is a **scalar loss**.
- The index $r$ is a reaction index, $j$ a component index, $i$ a cell index
  (overload resolution I1 of §2).
- $\sigma$ denotes the Stefan–Boltzmann constant where it appears in physics,
  not a standard deviation; variances in the global-sensitivity section are
  written $\mathbb{V}[\cdot]$ to avoid the clash.

---

## 17.1 Motivation: why sensitivities

A pyrolysis simulation maps a finite set of material and boundary parameters
$\boldsymbol{\theta}$ — Arrhenius pre-exponentials $A_r$ and activation energies
$E_r$, heats of reaction $h_r$, component densities $\rho_j$, conductivities
$k_j$, emissivities $\varepsilon_j$, incident-flux histories, convective
coefficients $h_{\text{conv}}$ — to time-resolved fields and to derived scalar
observables: surface temperature $T_s(t)$, mass-loss rate (MLR), peak heat
release rate, time to ignition, residual char fraction. The **sensitivity** of
an output to a parameter,

$$
S_{kq} \;=\; \frac{\partial y_k}{\partial \theta_q},
$$

is the workhorse quantity for three distinct tasks.

**(1) Inverse parameter estimation (calibration).** Cone-calorimeter or TGA
data fix observables $\mathbf{y}^{\text{data}}$; one seeks $\boldsymbol{\theta}$
minimizing a misfit
$\ell(\boldsymbol{\theta})=\tfrac12\lVert\mathbf{y}(\boldsymbol{\theta})-\mathbf{y}^{\text{data}}\rVert_W^2$.
Gradient-based optimizers (Gauss–Newton, L-BFGS, Levenberg–Marquardt) require
$\nabla_{\boldsymbol{\theta}}\ell$, which is exactly an adjoint or a
Jacobian-vector contraction of $S_{kq}$. The shipped inverse-analysis workflow
(`examples/inverse_analysis/wp_objective.jl`, `wp_forward.jl`) fits White-Pine
kinetics and properties to multi-flux cone data this way.

**(2) Uncertainty quantification (UQ).** Given parameter uncertainties
(covariances from instrument calibration and per-curve regression), the
first-order propagated output covariance is
$\boldsymbol{\Sigma}_y \approx S\,\boldsymbol{\Sigma}_\theta\,S^{\mathsf T}$,
the *delta method*. Variance-based (Sobol) decomposition (§17.5) goes beyond
this linearization to apportion output variance to parameter groups.

**(3) Optimization and design.** Maximizing fire-retardant performance,
selecting char-forming additives, or designing a TGA heating program are all
gradient-amenable once $\nabla_{\boldsymbol{\theta}}\ell$ is available.

Pyrolysis.jl provides two complementary *local* derivative engines and a
*global* harness:

| Mode | Engine | Cost scaling | Best when |
|---|---|---|---|
| Forward (local) | `forward_sensitivity` / ForwardDiff | $\mathcal{O}(p)$ solves' worth | small $p$ ($\lesssim 20$), many outputs |
| Adjoint (local) | `adjoint_sensitivity` / Zygote + SciMLSensitivity | $\mathcal{O}(1)$ in $p$ per scalar loss | large $p$, scalar loss (inverse problems) |
| Global (variance) | parameter-sweep + surrogate (Sobol) | $\mathcal{O}(N_{\text{samples}})$ solves | apportioning variance, ranking groups |

Both local engines are *wrapper-based*: they differentiate the existing
`solve(problem)` end to end rather than rewriting the residual. This is the
single most important design decision in this chapter and is treated in §17.4.

---

## 17.2 Forward-mode local sensitivity (ForwardDiff)

### 17.2.1 The differentiated map

`forward_sensitivity` computes the full Jacobian of a vector output with respect
to the parameters by forward-mode algorithmic differentiation of the composite
map

$$
\boldsymbol{\theta}
\;\xmapsto{\;p_{\text{inject}}\;}\;
\text{problem}'
\;\xmapsto{\;\texttt{solve}\;}\;
\text{solution}
\;\xmapsto{\;\mathcal{O}\;}\;
\mathbf{y}.
$$

The extension wraps these three stages into a single closure and hands it to
`ForwardDiff.jacobian`:

```julia
wrapped = θ -> output_fn(solve(p_inject(base_problem, θ); progress=false, solve_kwargs...))
output   = wrapped(θ)
∂out_∂θ  = ForwardDiff.jacobian(wrapped, θ)
return (output, ∂out_∂θ)
```

Mathematically the returned matrix is

$$
\big[\partial\mathbf{y}/\partial\boldsymbol{\theta}\big]_{kq}
=\frac{\partial\, \mathcal{O}_k\big(\texttt{solve}(p_{\text{inject}}(\text{base},\boldsymbol{\theta}))\big)}{\partial \theta_q},
\qquad k=1\dots m,\;\; q=1\dots p ,
$$

evaluated at the supplied base point $\boldsymbol{\theta}$. The signature is
`forward_sensitivity(base_problem, θ, p_inject; output_fn=identity, solve_kwargs...)`;
`p_inject` is a positional closure and `output_fn` defaults to identity (return
the whole solution object, which ForwardDiff then differentiates element-wise).

> **Implementation note.** Forward-mode sensitivity is implemented in
> `Pyrolysis.forward_sensitivity`,
> `ext/PyrolysisSciMLSensitivityExt.jl` (lines 30–48), attached to the stub
> `forward_sensitivity` declared in `src/sensitivity_stubs.jl`. The stub raises a
> `MethodError`-style error until `using SciMLSensitivity` (and `Zygote`) loads
> the extension. ForwardDiff is the only differentiation engine used by the
> *forward* path; SciMLSensitivity is a load trigger, not an active component
> here.

### 17.2.2 Dual propagation end to end

Forward mode works because the package carries the scalar element type `Tu`
of the state through every cache and every kernel. `ForwardDiff.jacobian`
seeds $\boldsymbol{\theta}$ with $p$-component dual numbers
$\theta_q + \sum_{s} \delta_{qs}\,\epsilon_s$ (one partial $\epsilon_s$ per
parameter; in practice ForwardDiff chunks the $p$ partials). Each dual
satisfies the algebra
$f(\theta+\dot\theta\,\epsilon) = f(\theta) + f'(\theta)\,\dot\theta\,\epsilon$,
so propagating duals through the solve simultaneously evaluates the output and
its directional derivatives. For the chain to be exact, no stage may drop the
dual part or silently coerce to `Float64`. The package satisfies this:

- **State cache.** `StateCache{NC,Tu,Mode}` is parameterized on `Tu`; for AD the
  workspace is built with `Tu = ForwardDiff.Dual{Tag,Float64,N}` (§15,
  `build_ad_workspace`). The column cross-section area is stored as
  `A_section::Base.RefValue{Tu}` so that the $\chi$-driven area carries
  sensitivities. (See `src/Problem/state_cache.jl`.)
- **State-derived geometry.** In ALE mode, cell volumes
  $V_i=\lvert z_{i+1}-z_i\rvert\,A_{\text{section}}$, cell centers, and face
  positions are recomputed from the dual-typed node positions in
  `update_state_geometry!(cache, mesh, ::ALE)`; the mutable mesh is read-only
  (connectivity only). Geometry sensitivities therefore propagate (§11, The ALE
  Framework).
- **Property evaluation.** Property functions are callable on `Real`, hence on
  duals. The analytical temperature derivatives of §4.5 — e.g. Arrhenius
  $\mathrm{d}p/\mathrm{d}T = p(T)\,E/(R_g T^2)$ and the unrolled polynomial
  derivative (`derivative`, `src/Materials/property_functions.jl:283–300`) — are
  used by the *Jacobian* assembly, but for *forward sensitivity* ForwardDiff
  differentiates the property's own arithmetic, so even `TableProperty` and
  `CallableProperty` (which only have finite-difference $\mathrm{d}p/\mathrm{d}T$
  internally) are differentiated *exactly* with respect to $T$ through their
  interpolation arithmetic.
- **Boundary-temperature Newton solve.** The surface-temperature solve
  $f(T_s)=0$ (§12.2) detects dual inputs and returns a dual $T_s$; the implicit
  function theorem is realized automatically because Newton's iterate
  $T_s\leftarrow T_s - f/f'$ is itself differentiable arithmetic once converged.
  (See `solve_boundary_temperature`, `src/Discretization/boundary_temperature.jl`.)
- **Kinetics.** `reaction_rate` is written with
  `one(promote_type(eltype(ξ), typeof(T)))` and `_fast_pow` so that the
  concentration product, tanh gates, and depletion limiter all stay in the
  promoted dual type (`src/Physics/kinetics.jl:75–111`).

The time integrator (`KenCarp4` by default, with
`autodiff=AutoFiniteDiff()` for its *internal* Newton Jacobian and the
default sparse factorization for the linear solve) does **not** need to be
dual-aware: ForwardDiff differentiates the *output of the converged solve*, not
the integrator's internal Jacobian. OrdinaryDiffEq's `solve` is itself a
differentiable program, and SciMLSensitivity supplies the dual-compatible
solve path so that `ForwardDiff.Dual` parameters flow through the time loop.

### 17.2.3 What the derivative actually is

For an output $y=\mathcal{O}(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta})$
evaluated at the final state, the total derivative obeys the chain rule

$$
\frac{\mathrm{d}y}{\mathrm{d}\theta_q}
=\frac{\partial\mathcal{O}}{\partial\boldsymbol{\theta}}\bigg|_q
+\frac{\partial\mathcal{O}}{\partial\mathbf{u}}\cdot
\frac{\partial\mathbf{u}(t_{\text{end}})}{\partial\theta_q}.
$$

The second factor is the *state sensitivity*
$\mathbf{s}_q(t)\equiv\partial\mathbf{u}(t)/\partial\theta_q$, which (formally)
satisfies the **forward sensitivity ODE** obtained by differentiating the
residual $\dot{\mathbf{u}}=\mathbf{f}(\mathbf{u},\boldsymbol{\theta},t)$:

$$
\dot{\mathbf{s}}_q
=\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\,\mathbf{s}_q
+\frac{\partial\mathbf{f}}{\partial\theta_q},
\qquad
\mathbf{s}_q(0)=\frac{\partial\mathbf{u}_0}{\partial\theta_q}.
\tag{17.1}
$$

Here $\partial\mathbf{f}/\partial\mathbf{u}$ is precisely the structured
Jacobian of §15. Pyrolysis.jl does **not** integrate (17.1) explicitly; the
forward-mode wrapper instead lets ForwardDiff carry the dual partials through
the existing time loop, which is algebraically equivalent to solving (17.1) by
the *same* discretization the integrator uses for the primal — i.e. it inherits
the integrator's order and step-size control automatically, including the
contribution of $\partial\mathbf{u}_0/\partial\theta_q$ when a parameter
appears in the initial condition (e.g. an initial density). Equation (17.1) is
the conceptual object; the dual-propagation implementation is the realized one.

### 17.2.4 Cost and when to use forward mode

Forward mode evaluates $p$ directional derivatives in lockstep with the primal.
Its cost is roughly that of $1+p$ primal solves (ForwardDiff chunks reduce the
constant but not the linear-in-$p$ scaling). It is therefore the right choice
when:

- the parameter count is small ($p\lesssim 20$ — a handful of reactions and
  property coefficients), and
- many outputs are wanted at once (the *full* $m\times p$ Jacobian comes for
  free; reverse mode would need $m$ passes for $m$ scalar outputs).

For a typical single-reaction calibration with $\boldsymbol{\theta}=(A,E,h)$ and
outputs $\{T_s(t_k),\text{MLR}(t_k)\}$, forward mode delivers the entire
sensitivity matrix in three dual-augmented solves' worth of work — ideal for
Gauss–Newton, whose normal equations need $S$ itself, not just
$S^{\mathsf T}\mathbf{r}$.

---

## 17.3 Adjoint-mode local sensitivity (Zygote + SciMLSensitivity)

### 17.3.1 The differentiated map and the gradient

When the quantity of interest is a *scalar* loss $\ell(\boldsymbol{\theta})$
(the inverse-problem case), reverse mode computes the entire gradient
$\nabla_{\boldsymbol{\theta}}\ell$ at a cost essentially independent of $p$.
`adjoint_sensitivity` wraps the loss-producing solve and differentiates it with
Zygote:

```julia
wrapped = θ -> loss_fn(solve(p_inject(base_problem, θ);
                             progress=false, sensealg=sensealg, solve_kwargs...))
loss = wrapped(θ)
grad = first(Zygote.gradient(wrapped, θ))
return (loss, grad)
```

The returned gradient is
$\big[\nabla_{\boldsymbol{\theta}}\ell\big]_q=\partial\ell/\partial\theta_q$.
The signature is
`adjoint_sensitivity(base_problem, θ, p_inject, loss_fn; sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP()), solve_kwargs...)`.

> **Implementation note.** Adjoint sensitivity is
> `Pyrolysis.adjoint_sensitivity`, `ext/PyrolysisSciMLSensitivityExt.jl`
> (lines 53–74). The `sensealg` is **forwarded to `solve` as a first-class
> keyword** — `solve` accepts `sensealg` precisely so SciMLSensitivity can
> intercept the ODE solve and augment it with adjoint dynamics without tripping
> `solve`'s unknown-keyword guard (`src/Solver/solve.jl`). Zygote provides the
> reverse-mode pullback over the *non-ODE* parts of the chain (the
> `p_inject`/`loss_fn` glue); SciMLSensitivity provides the differentiable ODE
> solve in between.

### 17.3.2 The continuous adjoint and the augmented ODE

The mathematics SciMLSensitivity realizes is the *continuous adjoint*. Let the
loss be an integral (and/or endpoint) functional of the trajectory,

$$
\ell=\int_{0}^{t_{\text{end}}} g\big(\mathbf{u}(t),\boldsymbol{\theta},t\big)\,\mathrm{d}t
\;+\;\Phi\big(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta}\big),
$$

subject to $\dot{\mathbf{u}}=\mathbf{f}(\mathbf{u},\boldsymbol{\theta},t)$.
Introducing the adjoint (costate) $\boldsymbol{\lambda}(t)\in\mathbb{R}^{N_{\text{state}}}$,
stationarity of the Lagrangian gives the **backward adjoint ODE**

$$
\dot{\boldsymbol{\lambda}}
=-\left(\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\right)^{\!\mathsf T}\boldsymbol{\lambda}
-\left(\frac{\partial g}{\partial\mathbf{u}}\right)^{\!\mathsf T},
\qquad
\boldsymbol{\lambda}(t_{\text{end}})=\left(\frac{\partial\Phi}{\partial\mathbf{u}}\right)^{\!\mathsf T}\bigg|_{t_{\text{end}}},
\tag{17.2}
$$

integrated from $t_{\text{end}}$ back to $0$, after which the gradient is the
quadrature

$$
\frac{\mathrm{d}\ell}{\mathrm{d}\theta_q}
=\frac{\partial\Phi}{\partial\theta_q}
+\int_{0}^{t_{\text{end}}}
\left(
\boldsymbol{\lambda}^{\mathsf T}\frac{\partial\mathbf{f}}{\partial\theta_q}
+\frac{\partial g}{\partial\theta_q}
\right)\mathrm{d}t .
\tag{17.3}
$$

The key structural facts:

- (17.2) is a single ODE of the same dimension as the state, **independent of
  $p$**. One backward solve yields the full gradient via the $p$ cheap
  quadratures (17.3). This is why reverse mode dominates when $p$ is large.
- The operator $(\partial\mathbf{f}/\partial\mathbf{u})^{\mathsf T}$ is the
  transpose of the structured Jacobian of §15. SciMLSensitivity never forms it
  densely; with `autojacvec=ReverseDiffVJP()` it evaluates *vector–Jacobian
  products* $\boldsymbol{\lambda}^{\mathsf T}(\partial\mathbf{f}/\partial\mathbf{u})$
  by reverse-mode AD of the residual, matching the matrix-free spirit of the
  `MatrixFreeSolve` strategy.
- The default `InterpolatingAdjoint` stores (and interpolates) the forward
  trajectory so the backward pass can re-evaluate $\mathbf{f}$ and its VJPs at
  the correct states; this trades memory for the avoidance of a backward
  re-solve of the primal.

### 17.3.3 sensealg selection

The `sensealg` keyword chooses the adjoint algorithm. The default is
`InterpolatingAdjoint(autojacvec = ReverseDiffVJP())`. The practical menu (all
from SciMLSensitivity) and their trade-offs:

| `sensealg` | Memory | Robustness | Notes |
|---|---|---|---|
| `InterpolatingAdjoint` (default) | high (stores forward sol) | high | recompute-free backward pass |
| `QuadratureAdjoint` | high | high | separate quadrature for (17.3); fast for few params |
| `BacksolveAdjoint` | low | lower (can drift for stiff Arrhenius) | re-solves primal backward; risky here |
| `GaussAdjoint` | moderate | high | checkpointed quadrature |
| `ForwardDiffSensitivity` | — | exact | reverts to forward mode; small $p$ only |

For pyrolysis the reactions are stiff and the trajectory can vary by orders of
magnitude across a decomposition front, so `BacksolveAdjoint` is generally
discouraged (its backward re-solve can leave the stiff manifold); the
interpolating or Gauss variants are preferred. The `autojacvec` choice
(`ReverseDiffVJP`, `ZygoteVJP`, `EnzymeVJP`, or `TrackerVJP`) sets how the VJP
$\boldsymbol{\lambda}^{\mathsf T}\partial\mathbf{f}/\partial\mathbf{u}$ is taken;
`ReverseDiffVJP` is robust for the scalar-heavy residual kernels here.

### 17.3.4 Cost and when to use adjoint mode

Reverse mode costs $\mathcal{O}(1)$ scalar-loss gradients independent of $p$
(one forward solve + one backward adjoint solve), at the price of higher memory
and the need for a fully reverse-differentiable solve path. Use it when:

- the loss is a single scalar (sum-of-squares misfit, integrated MLR error),
  and
- the parameter count is moderate-to-large (many property coefficients across
  several reactions, distributed-activation surrogates, spatially-varying
  initial fields).

The dividing line between §17.2 and §17.3 is essentially $p$ versus $m$: forward
mode for "few parameters, many outputs," adjoint mode for "many parameters, one
loss." For the in-between regime the wrapper lets a user simply swap
`forward_sensitivity` ↔ `adjoint_sensitivity` with the *same* `p_inject`.

---

## 17.4 The `p_inject` wrapper pattern and parameter handles

### 17.4.1 Why wrap rather than rewrite

Both engines avoid touching the residual signature. The residual
$\texttt{residual!}(du,u,\text{def},t,ws)$ takes an *immutable* `ProblemDef` and
a *mutable* `Workspace` (§15.2); it has no explicit parameter argument. Rather
than thread a parameter vector $\boldsymbol{\theta}$ through `def` and every
operator, sensitivity is obtained by differentiating the *construction* of the
problem followed by the solve. The user supplies a closure

$$
p_{\text{inject}}:\ (\text{base\_problem},\boldsymbol{\theta})\ \longmapsto\ \text{problem}',
$$

that materializes a *fresh* `PyrolysisProblem` with $\boldsymbol{\theta}$ baked
into the relevant fields, and the engine differentiates
$\boldsymbol{\theta}\mapsto\texttt{solve}(p_{\text{inject}}(\dots))$.

This has three consequences worth stating precisely:

1. **No invasive refactor.** The residual, the operators, the Jacobian backends,
   and the state layout are untouched. The cost is that AD must re-run *problem
   construction* (cheap: building a `Material` and `BoundaryConditionSet`) on
   every evaluation; the solve dominates anyway.
2. **Any constructible parameter is differentiable.** Because $p_{\text{inject}}$
   builds the problem, *anything* that flows into a `Material`,
   `BoundaryConditionSet`, mesh, or initial condition is a valid handle — there
   is no fixed parameter registry.
3. **The base problem must not be mutated.** `p_inject` returns a *new*
   problem; mutating `base_problem` in place would break the purity that lets
   sweeps and AD share immutable state (the same contract the parameter-sweep
   harness enforces, §17.5).

### 17.4.2 Parameter handles

Common $\boldsymbol{\theta}$ entries and how `p_inject` injects them:

- **Kinetic parameters** $A_r$, $E_r$. Rebuild the reaction(s) with new
  pre-exponential/activation energy and assemble a new `Material`. The rate
  $r=A\,e^{-E/(R_gT)}\prod_j\xi_j^{n_{r,j}}$ is analytic in both $A$ (linear) and
  $E$ (through the exponential), so $\partial r/\partial A = r/A$ and
  $\partial r/\partial E = -r/(R_gT)$ are exact under AD. Because $A$ and $E$ are
  strongly correlated along the Arrhenius "compensation" direction, the
  sensitivity matrix is typically ill-conditioned in $(A,E)$ — see §17.7.
- **Heat of reaction** $h_r$ (or a `LinearProperty` $a+bT$, or a
  `StateDependentProperty` for moisture sorption). It enters only the energy
  source $Q_{\text{rxn}}=-\sum_r h_r r_r$ linearly, so
  $\partial Q_{\text{rxn}}/\partial h_r = -r_r$ exactly. With the sign
  convention $h>0$ endothermic, increasing $h_r$ *cools* the cell, lowering
  $T_s$ — the sensitivity $\partial T_s/\partial h_r<0$.
- **Densities / conductivities / emissivities** $\rho_j,k_j,\varepsilon_j$.
  These feed the mixture rules of §5 and the boundary energy balance; their
  derivatives propagate through the effective-property closures (analytical
  mixing derivatives for PARALLEL/SERIES/WEIGHTED, AD for BRUGGEMAN). A
  bulk-density change also shifts the initial dry-solid mass $m_{\text{dry},i}$
  (and hence $\chi$, §10), so density is a particularly broad-reaching handle.
- **Boundary-condition parameters** — incident flux level, convective
  coefficient $h_{\text{conv}}$, ambient $T_\infty$, view factor $F$,
  emissivity in $q_{\text{rad}}=F\varepsilon\sigma(T_\infty^4-T_s^4)$. These
  enter the surface energy balance solved for $T_s$; their sensitivities pass
  through the implicit Newton solve as in §17.2.2.

### 17.4.3 Worked structure of a forward-sensitivity call

The canonical usage that the extension differentiates is:

```julia
using Pyrolysis, SciMLSensitivity, Zygote   # loads the extension

θ0 = [A0, E0, h0]                            # base parameter point

function p_inject(prob, θ)
    rxn  = Reaction(:decomp, 1 => 2; A = θ[1], E = θ[2], h = θ[3])
    mat  = Material(prob.material; reactions = (rxn,))
    return PyrolysisProblem(prob; material = mat)
end

output_fn(sol) = sol.surface_T               # vector QoI over saved times

(y, S) = forward_sensitivity(base_problem, θ0, p_inject; output_fn = output_fn)
# y :: surface-T trace; S[k,q] = ∂(T_s at t_k)/∂θ_q
```

`S` is the $m\times 3$ sensitivity matrix used directly in the Gauss–Newton
normal equations $S^{\mathsf T}WS\,\Delta\boldsymbol{\theta}
= -S^{\mathsf T}W\,\mathbf{r}$ for the residual
$\mathbf{r}=\mathbf{y}-\mathbf{y}^{\text{data}}$.

---

## 17.5 Global / variance-based sensitivity (Sobol)

### 17.5.1 What global sensitivity adds

Local sensitivities $S_{kq}$ are first-order Taylor coefficients at one point.
When parameter uncertainties are large or the response is nonlinear, one wants
*variance-based* (global) indices that hold over the whole input distribution.
Treat the parameters as random $\boldsymbol{\Theta}$ with a joint distribution
and a scalar QoI $Y=\mathcal{Q}(\boldsymbol{\Theta})$. The Sobol/ANOVA
decomposition apportions the output variance $\mathbb{V}[Y]$ to subsets of
inputs. For an input (or group) $\theta_q$ the **first-order** and **total**
indices are

$$
S_q^{\,1}=\frac{\mathbb{V}_{\theta_q}\!\big(\mathbb{E}_{\sim q}[Y\mid\theta_q]\big)}{\mathbb{V}[Y]},
\qquad
S_q^{\,T}=1-\frac{\mathbb{V}_{\sim q}\!\big(\mathbb{E}_{\theta_q}[Y\mid\boldsymbol{\theta}_{\sim q}]\big)}{\mathbb{V}[Y]},
\tag{17.4}
$$

where $\sim q$ denotes "all parameters except $q$." $S_q^1$ measures the variance
explained by $\theta_q$ acting alone; $S_q^T$ adds all interactions involving
$\theta_q$. They bracket the importance: $S_q^1\le S_q^T$, with equality iff
$\theta_q$ has no interactions.

### 17.5.2 Design, sampling, and the surrogate

Pyrolysis.jl does not bundle a Sobol estimator into the core solver; instead the
global workflow is built on the **parameter-sweep harness** (§17.5.3) plus a
*surrogate* fitted to the swept solutions. The shipped reference implementation
is `examples/sensitivity_analysis/sensitivity_sobol.jl`, which for the
Douglas-fir cone problem:

1. Draws correlated parameter samples (16 uncertain parameters in 7
   *independent* groups — calibration biases, low/high-T virgin $c_p$, char
   $c_p$, heat of sorption, heats of reaction). Correlations are *within-group
   only*, so the standard Sobol formulas (17.4) apply at the group level without
   correlated-input extensions.
2. Transforms the correlated draws to independent standard normals (Cholesky
   inverse for Gaussian groups; Rosenblatt sequential conditioning for
   Student-$t$ groups).
3. Solves the pyrolysis problem at each sample (the expensive step — one full
   cone solve per sample) and reduces each solution to scalar QoIs.
4. Fits an adaptive-degree (2–4) **Polynomial Chaos Expansion (PCE)** surrogate
   to each QoI, selecting degree by leave-one-out $Q^2$ and using a sparse LASSO
   fit when the oversampling ratio is low.
5. Computes first-order and total *group* Sobol indices **analytically from the
   PCE coefficients** — because the PCE basis is orthonormal under the input
   measure, each Sobol index is a normalized sum of squared coefficients, so no
   second Monte-Carlo loop is needed.
6. Quantifies index uncertainty with a BCa bootstrap and assesses convergence on
   nested subsamples.

The PCE surrogate is the practical enabler: it converts an
$\mathcal{O}(N_{\text{samples}})$ *solve* budget into closed-form indices,
avoiding the $\mathcal{O}((p{+}2)N)$ solves of a raw Saltelli/Sobol pick-freeze
estimator. The degree-selection rationale is documented in
`examples/sensitivity_analysis/pce_degree_selection.md`; `compare_degrees.jl`
sweeps the degree explicitly.

### 17.5.3 Relation to the parameter-sweep harness

The variance design is executed through `parameter_sweep(problem, grid; modify, extract)`
(§ User Guide UG-12). Its contract dovetails with the sampling above:

$$
\texttt{modify}(\text{base\_problem},\,\mathbf{p}_s)\ \longmapsto\ \text{problem}_s,
\qquad
\texttt{extract}(\text{solution}_s)\ \longmapsto\ \mathbf{y}_s,
$$

for each sample point $\mathbf{p}_s$ in the design `grid`. `modify` plays the
same role as `p_inject` in §17.4 — it must build a *fresh* problem and must not
mutate the base. The harness dispatches the (embarrassingly parallel) solves
over `Distributed.pmap` when workers are available, giving each worker its own
mutable `Workspace` while sharing the immutable `ProblemDef`:

```julia
work = point -> extract(solve(modify(problem, point); progress=false, solve_kwargs...))
results = (nworkers() > 1) ? pmap(work, pts) : map(work, pts)
```

(`src/Solver/parameter_sweep.jl`.) The `extract` closure should return small
NamedTuples of scalar QoIs (peak $T_s$, integrated MLR, time-to-peak) to
minimize worker→head transfer; those scalars are exactly the $Y$ samples the PCE
surrogate consumes. Thus the *same* solve and the *same* immutable-problem
discipline underpin local sensitivity (§17.2–17.4), parameter sweeps, and global
Sobol analysis — one harness, three derivative products.

### 17.5.4 Local-vs-global consistency

The two views are linked by the *delta method*: near a base point, the
first-order propagated variance is
$\mathbb{V}[Y]\approx\sum_{q,q'}S_q S_{q'}\,\mathrm{Cov}(\theta_q,\theta_{q'})$
with $S_q=\partial Y/\partial\theta_q$ from §17.2/§17.3. When the response is
linear in the inputs over their uncertainty range, the Sobol $S_q^1$ equals the
normalized squared local sensitivity
$S_q^2\,\mathbb{V}[\theta_q]/\mathbb{V}[Y]$ and $S_q^T=S_q^1$ (no interactions).
Discrepancies between the local and global rankings are therefore *diagnostics
of nonlinearity/interaction*, not contradictions — a useful cross-check when
both are available.

---

## 17.6 Interaction with AD-clean design choices

Sensitivities are only meaningful if the map $\boldsymbol{\theta}\mapsto\mathbf{y}$
is (at least) continuously differentiable. Several modeling choices in
Pyrolysis.jl exist *specifically* to make this true, and they are precisely the
smooth replacements for the classical non-smooth pyrolysis closures. Each is
verified against the source below.

**(1) Smooth temperature gates.** The hard Heaviside window
$H(T-T_{\min})H(T_{\max}-T)$ would make $\partial r/\partial T$ a delta function
at the gate edges and destroy both the integrator's Newton Jacobian and any AD
derivative. The code replaces it with the product of tanh ramps

$$
r \;\mathrel{*}=\; \tfrac12\big(1+\tanh(T-T_{\min})\big)\,
                   \tfrac12\big(1+\tanh(T_{\max}-T)\big),
$$

applied only when the bounds are finite/positive (`src/Physics/kinetics.jl:102–108`).
This is $C^\infty$ in $T$, so all temperature sensitivities of the rate are
finite. The cost is a $\sim\!10\text{ K}$ leakage of reactivity outside the
nominal window (the ramp is not unit-step), discussed in §6.2; for sensitivity
purposes the benefit — globally finite $\partial r/\partial T$ — outweighs it.

**(2) Depletion limiter.** Near zero reactant the bare Arrhenius rate is stiff
and its parameter sensitivities blow up. Each reactant multiplies the rate by
the smooth factor

$$
\tanh(\xi/\xi_{\text{thr}}),
$$

with default threshold $\xi_{\text{thr}}=1.0\ \text{kg·m}^{-3}$
(`src/Physics/kinetics.jl:89–95`; reactant clamped via $\max(\xi,0)$ first).
The rate vanishes at depletion with a bounded $\partial r/\partial\xi$ for
every reaction order $n \ge 0$, keeping the reaction-rate Jacobian — and hence
the parameter sensitivities of MLR near burnout — bounded and well-defined.
The factor is $C^\infty$ in $\xi > 0$ with no threshold branch; the only
non-smooth point left is the $\max(\xi,0)$ clip exactly at $\xi=0$ (§17.7).

**(3) Midpoint-rule sensible enthalpy.** The sensible-enthalpy increment
$\Delta h(T_{\text{hi}},T_{\text{lo}})=c_p\!\big(\tfrac{T_{\text{hi}}+T_{\text{lo}}}2\big)(T_{\text{hi}}-T_{\text{lo}})$
(`src/Materials/components.jl:515–518`) uses a *single* $c_p$ evaluation at the
midpoint. It is smooth in both temperatures and in any property parameters
inside $c_p$, and — crucially — AD returns the derivative of *what the midpoint
rule actually computes*, not of an idealized exact integral. This consistency
means energy-balance-derived QoIs differentiate cleanly and match the
single-point cell-centered $c_p$ used elsewhere.

**(4) Smooth flux limiters and lateral-shrinkage derivative.** The TVD limiters
(Van Leer is the ALE default) are continuous, and the $\chi$-driven radial
strain rate uses an *automatic-differentiation* derivative of the
lateral-shrinkage law,
$\theta_A=\dfrac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\dfrac{\mathrm{d}\bar\chi}{\mathrm{d}t}$
with `dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)` and a guard
`A_ratio > eps()` (`src/Residual/residual.jl:559–563`). Because $\mathrm{law}'$
is taken by AD, *any differentiable* law contributes correct geometry
sensitivities — but a non-smooth law breaks this (see §17.7).

**(5) Harmonic-mean face averaging and tanh-based smoothing elsewhere.** The
distance-weighted harmonic mean used for face $k_f,\lambda_f,\kappa_f$ (§5.12)
and the smooth (rather than `min`/`max`) constructions in the gates keep the
discrete operators differentiable. The few genuinely non-smooth operations that
remain — `clamp` on porosity/temperature, `max` in the depletion clamp, and the
state rewrite at park events — are localized and are exactly the pitfalls
catalogued next.

> **Implementation note.** The smoothness-critical kernels are
> `reaction_rate` (`src/Physics/kinetics.jl`), `Δh` /
> `sensible_enthalpy` (`src/Materials/components.jl`), the limiter library
> (`src/Discretization/flux_limiters.jl`), and `_compute_theta_A`
> (`src/Residual/residual.jl`). Each was written with `tanh`/midpoint/AD
> derivatives *in place of* the textbook step/cutoff/finite-difference form so
> that the chain $\boldsymbol{\theta}\mapsto\mathbf{y}$ is differentiable.

---

## 17.7 Pitfalls

Even with an AD-clean residual, sensitivity computations can return wrong or
meaningless numbers if certain solver features are active or the parameter point
sits on a non-smooth manifold. The pitfalls below are specific to Pyrolysis.jl.

### 17.7.1 Warm-start artifacts

`solve(...; warm_start=prev_solution)` reuses the prior solve's *final* state as
the new initial condition (`_apply_warm_start!`, `src/Solver/solve.jl`). Inside
a sensitivity loop this introduces a **hidden parameter dependence that AD does
not see**: the warm-start state was produced at a *different* $\boldsymbol{\theta}$,
so $\partial\mathbf{u}_0/\partial\theta_q$ is silently set to zero even though
the true initial condition (had the perturbed problem been integrated from the
real start) would depend on $\boldsymbol{\theta}$. Differentiating a
warm-started solve therefore yields a derivative of the *wrong map*. Open
question 3 in the solver KB flags this explicitly. **Rule:** do not warm-start
solves that are inside `forward_sensitivity`/`adjoint_sensitivity` or a
sensitivity sweep; integrate each sample from its true initial condition.

### 17.7.2 Non-smooth lateral-shrinkage laws

The radial strain rate uses $\mathrm{law}'(\bar\chi)$ via AD (§17.6.4). Identity
(`NoChi`) and the smooth linear laws $1\mp k\bar\chi$ are differentiable, but a
piecewise or table-defined law with kinks makes $\mathrm{law}'$ discontinuous at
the breakpoints. A sensitivity evaluated *at* a breakpoint is then ill-defined
(left/right derivatives differ), and AD returns whichever branch the
floating-point comparison happens to take. The guard `A_ratio > eps()` further
*zeroes* $\theta_A$ — and thus its sensitivity — when the cross-section nearly
vanishes, a hard switch that is correct physically but non-differentiable.
**Rule:** use $C^1$ shrinkage laws for sensitivity studies, and avoid evaluating
at $\bar\chi$ values where the law has a corner or where $\mathrm{law}(\bar\chi)\to 0$.

### 17.7.3 Park-event discontinuities

The depletion path (§14.3) parks a depleted cell: it merges the cell into a
neighbour and moves the freed slot above the live surface, rewriting the state
vector in place. The state keeps its size and the structured Jacobian keeps
working, but the rewrite is a **discrete event**: a perturbation of
$\boldsymbol{\theta}$ that shifts a park *time* produces a jump in the
trajectory, so $\mathbf{u}(t)$ is only piecewise smooth in the parameters. Forward
sensitivities are valid on horizons without park events; adjoint sensitivity,
which needs a differentiable solve path, cannot be run through a park event.

### 17.7.4 Arrhenius $(A,E)$ ill-conditioning

The pre-exponential and activation energy are nearly redundant along the
compensation direction $\delta\ln A \approx \delta E/(R_g T_{\text{peak}})$: a
joint perturbation that holds the rate fixed at the dominant reaction
temperature leaves the output almost unchanged. The sensitivity matrix $S$ is
consequently near-rank-deficient in the $(A,E)$ block, so Gauss–Newton normal
equations $S^{\mathsf T}WS$ are ill-conditioned. This is a property of the
*physics*, not a bug; mitigate by reparameterizing (e.g. fit $E$ and the rate at
a reference temperature instead of $A$ directly), by regularization
(Tikhonov/priors), or by using observables at *multiple* heating rates to break
the degeneracy.

### 17.7.5 Integrator-tolerance noise in finite-difference checks

When validating AD sensitivities against finite differences, remember the primal
solve is only accurate to the integrator tolerances (`abstol=1e-8`,
`reltol=1e-6` by default). A finite-difference step $\delta\theta$ smaller than
the solve noise yields garbage; a step too large incurs truncation error. AD
(forward or adjoint) avoids this entirely — the dual/adjoint derivative is
*exact* for the discretized solve at the chosen tolerances — which is precisely
why the package differentiates the solve rather than finite-differencing it. If a
finite-difference cross-check is needed, tighten `abstol`/`reltol` first.

### 17.7.6 Non-differentiable guards in the kinetics and clamps

The reactant clamp $\max(\xi,0)$, the temperature clamp to $[200,3000]\ \text{K}$
in the boundary Newton solve, and the porosity `clamp(\cdot,0,1)` are all
non-smooth at their boundaries. In normal operation states stay in the interior
and these are inactive, but a sensitivity evaluated exactly *at* a clamp boundary
(e.g. a fully depleted reactant, $\xi=0$) inherits a one-sided or zero
derivative. This is rarely an issue for QoIs integrated over time but can surface
in pointwise late-time MLR sensitivities near burnout.

---

## 17.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

Pyrolysis.jl's sensitivity capability is structurally different from its
condensed-phase peers, which historically expose parameter sensitivity only
*externally* (through driver scripts and finite differences) rather than through
differentiable internals.

**ThermaKin2Ds.** ThermaKin couples a forward pyrolysis solver to inverse
analysis via external optimization (e.g. Stochastic Hill-Climbing / Shuffled
Complex Evolution over the kinetic and property parameters), evaluating the
forward model as a black box and forming sensitivities by repeated forward runs
(finite differences). Mapping notation: ThermaKin's stoichiometric coefficients
$\theta_i^j$ are our yields $\nu_{r,j}$, and its heats of reaction carry the
*opposite* sign ($h>0$ exothermic) which our $Q_{\text{rxn}}=-\sum_r h_r r_r$
converts to (overload H1, §2). Pyrolysis.jl's advantage is that the *same*
parameter handles ($A_r,E_r,h_r,\rho_j,k_j$) are differentiated *analytically
through the solve* by AD, so gradients are exact at the integrator tolerance and
cost $\mathcal{O}(p)$ (forward) or $\mathcal{O}(1)$ (adjoint) rather than $p{+}1$
black-box forward runs per gradient.

**Gpyro.** Gpyro's calibration (Lautenberger & Fernandez-Pello) is built on a
*genetic algorithm* (GA) that minimizes the mismatch between simulated and
measured MLR/temperature over a population of parameter sets — derivative-free by
design. Gpyro maps the pre-exponential $Z$ → our $A_r$ and the conversion-based
nth-order kinetics to our mass-concentration form (§6.9). A GA explores a global
landscape without gradients (robust to the $(A,E)$ degeneracy and to non-smooth
responses) but needs many forward solves and gives no local sensitivity matrix.
Pyrolysis.jl is complementary: AD gradients drive fast *local* refinement
(Gauss–Newton/L-BFGS) once a GA or sweep has found a basin, and the
variance-based Sobol harness (§17.5) supplies the global parameter ranking that
a GA does not.

**FDS (solid phase).** The FDS solid-phase pyrolysis model recovers kinetic
parameters from a reference TGA temperature and peak rate (the
$A,E$-from-$T_{\text{peak}}$ recipe, §6.9) rather than offering a built-in
sensitivity/adjoint capability; FDS's broader uncertainty treatment is the
verification-and-validation accounting of its guides, not parametric
differentiation of the solid solver. Pyrolysis.jl differs in providing
first-class, AD-exact local sensitivities and a surrogate-based global
decomposition as part of the solver API.

In short: where the reference codes treat the forward model as a black box and
obtain sensitivities by repeated evaluation (finite difference) or avoid them
entirely (GA), Pyrolysis.jl exposes the forward map as a differentiable program
and offers both exact local gradients (forward and adjoint) and variance-based
global indices, all sharing one immutable-problem / mutable-workspace harness.

---

## 17.9 Limitations and open issues

Drawing on the open questions recorded in the solver and kinetics knowledge
bases:

1. **Adjoint sensitivity across park events (open).** A park is a discrete
   rewrite of the state, so the trajectory is only piecewise smooth across
   events and adjoint sensitivity cannot run through a depleting
   (`handle_depletion=true`) horizon that contains one. Reconciling the event
   with a differentiable solve path is an open design problem (Solver KB Q2).

2. **Warm-start in sensitivity loops (open).** Whether — and how — reusing a
   prior solution's final state as a fresh initial condition perturbs sensitivity
   estimates is unquantified (Solver KB Q3). The current guidance is to avoid
   warm-starting inside differentiated solves (§17.7.1).

3. **No native forward-sensitivity ODE / adjoint.** The wrapper approach
   re-runs problem construction and relies on SciMLSensitivity's generic
   machinery rather than a hand-coded sensitivity/adjoint exploiting the
   structured Jacobian and its $\mathcal{O}(N)$ couplings (§15.9). A native
   adjoint that transposes the structured operators directly would be faster and
   lower-memory; the API surface notes this as a planned future replacement.

4. **Non-smoothness residue.** A handful of guards (`max`/`clamp`, the
   `A_ratio > eps()` switch, hard temperature/
   porosity clamps) remain non-differentiable at their boundaries (§17.7.2,
   §17.7.6). They are inactive in the interior of normal operation but can
   contaminate sensitivities evaluated exactly on the boundary; a fully smooth
   reformulation of these guards is not yet in place.

5. **Parameter handles are unregistered.** Because injection is by user closure
   (`p_inject`/`modify`), there is no enforced catalogue of "supported"
   parameters or automatic correctness check that the injected field actually
   flows differentiably into the residual; responsibility for a clean,
   non-mutating, $C^1$ injection rests with the user (§17.4.1).

6. **Global-sensitivity tooling lives in examples, not the core.** Sobol/PCE is
   implemented in `examples/sensitivity_analysis/`, not as a packaged function;
   the design and degree-selection choices (PCE order, BCa bootstrap,
   group-correlation handling) are reference recipes rather than a stable API.

Cross-references: the structured Jacobian whose transpose underlies the adjoint
is §15 (Temporal Integration & Structured Jacobian); the smooth kinetics gates
and depletion limiter are §6 (Reaction Kinetics); the mixture-property
derivatives are §5 (Effective Properties); the lateral-shrinkage law and $\chi$
dynamics are §10 (Volume Change & Lateral Shrinkage); the ALE geometry whose
dual-typed recompute makes geometry sensitivities possible is §11 (The ALE
Framework); and the conservation diagnostics used to sanity-check
sensitivity-driven calibrations are §16 (Conservation Diagnostics). Validation
of sensitivity-driven fits against experiment and reference codes is §18
(Verification & Validation).
