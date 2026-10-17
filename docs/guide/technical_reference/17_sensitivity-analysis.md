# 17. Sensitivity Analysis (Theory)

This chapter develops the mathematics of parametric sensitivity for the
Pyrolysis.jl solver: how a derivative of a simulation output with respect to
model parameters is obtained, how accurate it is, and how it supports inverse
parameter estimation, uncertainty quantification, and optimization. We treat
both *local* sensitivity — central finite differences of the entire `solve`
chain — and *global*, variance-based sensitivity (Sobol indices computed through
the parameter-sweep harness). The solve path is `Float64` throughout: the
material constructors convert every parameter to `Float64`, and the integrator's
right-hand side, Jacobian and solution containers are `Float64`-typed. Neither a
`ForwardDiff.Dual` nor a reverse-mode pullback can therefore pass from a
parameter to an output, and local derivatives are finite differences of whole
solves (§17.3 explains what exact derivatives would require). Neither flavor has
a dedicated API; both are short user scripts around `solve`, and the
mathematics below is what those scripts compute.

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
$\nabla_{\boldsymbol{\theta}}\ell$ or the sensitivity matrix $S$ itself. The
shipped inverse-analysis workflow
(`examples/inverse_analysis/inverse_white_pine.jl`, with `wp_objective.jl` and
`wp_forward.jl`) minimizes such a misfit for White-Pine kinetics and properties
against multi-flux cone data, derivative-free (`Fminbox(NelderMead())`), with a
`FiniteDiff` Hessian at the optimum for the standard errors.

**(2) Uncertainty quantification (UQ).** Given parameter uncertainties
(covariances from instrument calibration and per-curve regression), the
first-order propagated output covariance is
$\boldsymbol{\Sigma}_y \approx S\,\boldsymbol{\Sigma}_\theta\,S^{\mathsf T}$,
the *delta method*. Variance-based (Sobol) decomposition (§17.5) goes beyond
this linearization to apportion output variance to parameter groups.

**(3) Optimization and design.** Maximizing fire-retardant performance,
selecting char-forming additives, or designing a TGA heating program are all
gradient-amenable once $\nabla_{\boldsymbol{\theta}}\ell$ is available.

Pyrolysis.jl supports one *local* derivative method and a *global* harness:

| Mode | Method | Cost | Best when |
|---|---|---|---|
| Local | central differences of `solve` (§17.2) | $2p+1$ solves | small $p$ (a handful of parameters), any number of outputs |
| Global (variance) | parameter sweep + surrogate (Sobol, §17.5) | $\mathcal{O}(N_{\text{samples}})$ solves | apportioning variance, ranking groups |

Both are *wrapper-based*: they evaluate the existing `solve(problem)` at chosen
parameter points rather than rewriting the residual. The closure that turns a
parameter vector into a problem is treated in §17.4.

---

## 17.2 Local sensitivity by central differences

### 17.2.1 The differentiated map

The local derivative is that of the composite map

$$
\boldsymbol{\theta}
\;\xmapsto{\;p_{\text{inject}}\;}\;
\text{problem}'
\;\xmapsto{\;\texttt{solve}\;}\;
\text{solution}
\;\xmapsto{\;\mathcal{O}\;}\;
\mathbf{y},
$$

approximated column by column with central differences,

$$
\big[\partial\mathbf{y}/\partial\boldsymbol{\theta}\big]_{:,q}
\;\approx\;
\frac{\mathbf{y}(\boldsymbol{\theta}+h_q\mathbf{e}_q)-\mathbf{y}(\boldsymbol{\theta}-h_q\mathbf{e}_q)}{2h_q},
\qquad q=1\dots p ,
$$

with $h_q=\max(\eta\,\lvert\theta_q\rvert,\,\eta_{\text{abs}})$ for a relative
step $\eta$ and an absolute floor $\eta_{\text{abs}}$ (FiniteDiff's `relstep` and
`absstep`; the floor defaults to $\eta$).
In practice `FiniteDiff.finite_difference_jacobian(y, θ, Val(:central); relstep = η)`
does the bookkeeping, with `y(θ) = O(solve(p_inject(base, θ)))` (User Guide
§13.3). Because the step is relative, parameters spanning many decades are
passed in a scaled form — $\log_{10}A_r$ rather than $A_r$ — and transformed
back inside $p_{\text{inject}}$.

### 17.2.2 Error balance and the noise floor

Write a computed output as $\tilde y(\theta)=y(\theta)+\epsilon(\theta)$, where
$y$ is the smooth output of the discretized model and $\epsilon$ the error the
time integration leaves in it. The central difference then carries two errors,

$$
\frac{\tilde y(\theta+h)-\tilde y(\theta-h)}{2h}-y'(\theta)
\;=\;\underbrace{\tfrac16\,y'''(\theta)\,h^{2}}_{\text{truncation}}
\;+\;\underbrace{\frac{\epsilon(\theta+h)-\epsilon(\theta-h)}{2h}}_{\text{noise}}
+\mathcal{O}(h^4).
$$

If $\epsilon$ were a smooth function of $\theta$ the noise term would be a
derivative too, and small. In general it is only piecewise smooth: the adaptive
step sequence, Newton iteration counts and Jacobian reuse decisions change
discretely as $\theta$ moves. A perturbation small enough to change none of
those decisions leaves $\epsilon$ smooth; one that changes a decision makes the
noise term of order $\lvert\epsilon\rvert/h$. Balancing that bound against
truncation gives the classical optimum
$h^\ast\sim\big(3\lvert\epsilon\rvert/\lvert y'''\rvert\big)^{1/3}$.

On the charring-slab example of User Guide §13.2 (40 cells, 50 kW m⁻²,
`reltol = 1e-7`) the estimates are stable over several decades of step: the
derivatives of the peak and 150-s surface temperatures, the mean surface
temperature and the final mass at $\eta=10^{-5}$–$10^{-7}$ agree with those at
$\eta=10^{-4}$ to $10^{-4}$ relative or better. Truncation error appears at
larger steps — up to 0.4 % at $\eta=10^{-3}$ and 29 % at $\eta=10^{-2}$, largest
for the 150-s temperature, which moves with the reaction front — and the error
against a tight reference falls with the tolerance, from $4\times10^{-3}$ K at
`reltol = 1e-5` to $2\times10^{-5}$ K at `1e-9`.

The same example with the condensed components at their default $\lambda=0$
behaves very differently (§17.7.7): 38 failed Newton iterations per solve
instead of none, a 0.35 K change in $T_s$ from a $10^{-6}$ relative change in
$A$, an error against the reference that stays at 0.2–0.4 K for every `reltol`
from $10^{-5}$ to $10^{-9}$, and derivatives that change sign from one step to
the next.

The practical rule: start from $\eta=10^{-4}$, confirm with $10\eta$ and
$\eta/10$, and treat disagreement as a sign that $y$ is not smooth in
$\boldsymbol{\theta}$ at this point (§17.7).

### 17.2.3 The state sensitivity being approximated

For an output $y=\mathcal{O}(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta})$
evaluated at the final state, the total derivative obeys the chain rule

$$
\frac{\mathrm{d}y}{\mathrm{d}\theta_q}
=\frac{\partial\mathcal{O}}{\partial\boldsymbol{\theta}}\bigg|_q
+\frac{\partial\mathcal{O}}{\partial\mathbf{u}}\cdot
\frac{\partial\mathbf{u}(t_{\text{end}})}{\partial\theta_q}.
$$

The second factor is the *state sensitivity*
$\mathbf{s}_q(t)\equiv\partial\mathbf{u}(t)/\partial\theta_q$, which satisfies
the **forward sensitivity ODE** obtained by differentiating the residual
$\dot{\mathbf{u}}=\mathbf{f}(\mathbf{u},\boldsymbol{\theta},t)$:

$$
\dot{\mathbf{s}}_q
=\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\,\mathbf{s}_q
+\frac{\partial\mathbf{f}}{\partial\theta_q},
\qquad
\mathbf{s}_q(0)=\frac{\partial\mathbf{u}_0}{\partial\theta_q}.
\tag{17.1}
$$

Here $\partial\mathbf{f}/\partial\mathbf{u}$ is precisely the structured
Jacobian of §15. Pyrolysis.jl does **not** integrate (17.1); the central
difference approximates $\mathbf{s}_q$ by the difference of two primal
trajectories, including the contribution of
$\partial\mathbf{u}_0/\partial\theta_q$ when a parameter appears in the initial
condition (e.g. an initial density), because $p_{\text{inject}}$ rebuilds the
initial condition too. Equation (17.1) is the object a native implementation
would integrate (§17.3.2).

### 17.2.4 Cost and parallelism

A central-difference Jacobian costs $2p$ perturbed solves plus the base solve,
independent of the number of outputs $m$: the full $m\times p$ matrix $S$ comes
from the same solves. That suits Gauss–Newton, whose normal equations need $S$
itself, and a handful of parameters. The $2p$ solves are independent;
`parameter_sweep` (§17.5.3) runs a list of perturbed parameter vectors with
$p_{\text{inject}}$ as its `modify` and distributes them over `Distributed`
workers (User Guide §13.5). A gradient of a scalar loss costs the same $2p$
solves, so for large $p$ a derivative-free or global method is usually the
better investment.

---

## 17.3 Why the solve is not differentiated automatically

### 17.3.1 What stops dual numbers and adjoints

The residual kernels are scalar-type-generic — they run on `ForwardDiff.Dual`
state for the Jacobian (§15) — but the path from a *parameter* to an output is
not:

1. **Construction.** Scalar properties go through
   `to_property(x::Real) = ConstantProperty(Float64(x))`
   (`src/Materials/property_functions.jl`), and `Reaction` stores $A$ and $E$ as
   `Float64` (`src/Materials/reactions.jl`). A dual-valued parameter fails with
   `MethodError: no method matching Float64(::ForwardDiff.Dual{…})` before any
   solve starts.
2. **Workspace and integrator.** `build_workspace` allocates every cache with
   `Tu = Float64`; the type-erased right-hand side (`ErasedRHS`) and Jacobian
   (`ErasedJacobian`, which accepts only an `AbstractMatrix{Float64}`) are
   `Float64`-typed so that the integrator compiles once for every problem (§15);
   `PyrolysisSolution` stores its traces as `Float64` arrays.
3. **No parameter vector.** The `ODEProblem` built by `solve` carries no
   parameter vector (`p = NullParameters()`); $\boldsymbol{\theta}$ enters only
   through closures over the constructed problem, where an adjoint
   sensitivity algorithm cannot see it.

### 17.3.2 What exact derivatives would require

Exact forward derivatives would integrate (17.1) alongside the primal, reusing
the structured Jacobian $\partial\mathbf{f}/\partial\mathbf{u}$ of §15 and
evaluating $\partial\mathbf{f}/\partial\theta_q$ by forward-mode AD of the
residual with respect to $\boldsymbol{\theta}$ alone. For a scalar loss
$\ell=\int_0^{t_{\text{end}}}g(\mathbf{u},\boldsymbol{\theta},t)\,\mathrm{d}t+\Phi(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta})$
and many parameters, the continuous adjoint is cheaper: the costate
$\boldsymbol{\lambda}(t)$ solves the backward ODE

$$
\dot{\boldsymbol{\lambda}}
=-\left(\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\right)^{\!\mathsf T}\boldsymbol{\lambda}
-\left(\frac{\partial g}{\partial\mathbf{u}}\right)^{\!\mathsf T},
\qquad
\boldsymbol{\lambda}(t_{\text{end}})=\left(\frac{\partial\Phi}{\partial\mathbf{u}}\right)^{\!\mathsf T}\bigg|_{t_{\text{end}}},
\tag{17.2}
$$

after which the gradient is the quadrature

$$
\frac{\mathrm{d}\ell}{\mathrm{d}\theta_q}
=\frac{\partial\Phi}{\partial\theta_q}
+\int_{0}^{t_{\text{end}}}
\left(
\boldsymbol{\lambda}^{\mathsf T}\frac{\partial\mathbf{f}}{\partial\theta_q}
+\frac{\partial g}{\partial\theta_q}
\right)\mathrm{d}t ,
\tag{17.3}
$$

one backward solve of the state's dimension for all $p$ parameters, with the
transpose of the structured Jacobian as its operator. Either route needs the
parameters lifted into a parameter vector the residual reads (item 3 above) and
a residual that is generic in the parameter element type (item 1). Both also
inherit the event caveats of §17.7.3: a park is a discrete state rewrite with
no sensitivity jump condition.

---

## 17.4 The `p_inject` wrapper pattern and parameter handles

### 17.4.1 Why wrap rather than rewrite

The residual $\texttt{residual!}(du,u,\text{def},t,ws)$ takes an *immutable*
`ProblemDef` and a *mutable* `Workspace` (§15.2); it has no explicit parameter
argument. Rather than thread a parameter vector $\boldsymbol{\theta}$ through
`def` and every operator, sensitivities are obtained by re-running *problem
construction* followed by the solve at each parameter point. The user supplies a
closure

$$
p_{\text{inject}}:\ (\text{base\_problem},\boldsymbol{\theta})\ \longmapsto\ \text{problem}',
$$

that materializes a *fresh* `PyrolysisProblem` with $\boldsymbol{\theta}$ baked
into the relevant fields, and the finite difference (or the Sobol sampler)
evaluates $\boldsymbol{\theta}\mapsto\texttt{solve}(p_{\text{inject}}(\dots))$.

This has three consequences worth stating precisely:

1. **No invasive refactor.** The residual, the operators, the Jacobian backends,
   and the state layout are untouched. Problem construction (building a
   `Material` and `BoundaryConditionSet`) is re-run at every evaluation; it is
   cheap next to the solve.
2. **Any constructible parameter can be perturbed.** Because $p_{\text{inject}}$
   builds the problem, *anything* that flows into a `Material`,
   `BoundaryConditionSet`, mesh, or initial condition is a valid handle — there
   is no fixed parameter registry.
3. **The base problem must not be mutated.** `p_inject` returns a *new*
   problem; mutating `base_problem` in place would make one evaluation leak into
   the next (the same contract the parameter-sweep harness enforces, §17.5).

### 17.4.2 Parameter handles

Common $\boldsymbol{\theta}$ entries and how they enter the model:

- **Kinetic parameters** $A_r$, $E_r$. Rebuild the reaction(s) with new
  pre-exponential/activation energy and assemble a new `Material`. The rate
  $r=A\,e^{-E/(R_gT)}\prod_j\xi_j^{n_{r,j}}$ is analytic in both $A$ (linear) and
  $E$ (through the exponential), with $\partial r/\partial A = r/A$ and
  $\partial r/\partial E = -r/(R_gT)$. Because $A$ and $E$ are strongly
  correlated along the Arrhenius "compensation" direction, the sensitivity
  matrix is typically ill-conditioned in $(A,E)$ — see §17.7.4. Perturb
  $\log_{10}A_r$, not $A_r$ (§17.2.1).
- **Heat of reaction** $h_r$ (or a `LinearProperty` $a+bT$, or a
  `StateDependentProperty` for moisture sorption). It enters only the energy
  source $Q_{\text{rxn}}=-\sum_r h_r r_r$ linearly, so
  $\partial Q_{\text{rxn}}/\partial h_r = -r_r$. With the sign convention
  $h>0$ endothermic, increasing $h_r$ *cools* the cell, lowering $T_s$ — the
  sensitivity $\partial T_s/\partial h_r<0$.
- **Densities / conductivities / emissivities** $\rho_j,k_j,\varepsilon_j$.
  These feed the mixture rules of §5 and the boundary energy balance. A
  bulk-density change also shifts the initial dry-solid mass $m_{\text{dry},i}$
  (and hence $\chi$, §10), so density is a particularly broad-reaching handle.
- **Boundary-condition parameters** — incident flux level, convective
  coefficient $h_{\text{conv}}$, ambient $T_\infty$, view factor $F$,
  emissivity in $q_{\text{rad}}=F\varepsilon\sigma(T_\infty^4-T_s^4)$. These
  enter the surface energy balance solved for $T_s$ (§12.2).

### 17.4.3 Worked structure of a local-sensitivity evaluation

```julia
using Pyrolysis, FiniteDiff

θ0 = [log10(A0), E0, h0]                     # base parameter point

function p_inject(prob, θ)
    rxn  = Reaction(:decomp, 1 => 2; A = 10.0^θ[1], E = θ[2], h = θ[3])
    m    = prob.material
    mat  = Material(; name = m.name, components = m.components, reactions = (rxn,),
                    mixing = m.mixing, lateral_shrinkage_law = m.lateral_shrinkage_law,
                    depletion_limiter = m.depletion_limiter, gas_viscosity = m.gas_viscosity)
    return PyrolysisProblem(; mesh = deepcopy(prob.mesh), material = mat,
                            bc_set = prob.bc_set, T_initial = prob.T_initial,
                            ξ_initial = prob.ξ_initial, tspan = prob.tspan,
                            dt_initial = prob.dt_initial, output_times = prob.output_times)
end

y(θ) = [solve(p_inject(base_problem, θ); progress = false).extras.total_mass[end]]

S = FiniteDiff.finite_difference_jacobian(y, θ0, Val(:central); relstep = 1e-4)
# S[k, q] = ∂y_k/∂θ_q
```

With a residual vector $\mathbf{r}=\mathbf{y}-\mathbf{y}^{\text{data}}$ in
place of $\mathbf{y}$, `S` is the $m\times p$ sensitivity matrix of the
Gauss–Newton normal equations $S^{\mathsf T}WS\,\Delta\boldsymbol{\theta}
= -S^{\mathsf T}W\,\mathbf{r}$.

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
Sobol analysis — and the same harness runs the perturbed solves of a
finite-difference Jacobian in parallel (§17.2.4).

### 17.5.4 Local-vs-global consistency

The two views are linked by the *delta method*: near a base point, the
first-order propagated variance is
$\mathbb{V}[Y]\approx\sum_{q,q'}S_q S_{q'}\,\mathrm{Cov}(\theta_q,\theta_{q'})$
with $S_q=\partial Y/\partial\theta_q$ from §17.2. When the response is
linear in the inputs over their uncertainty range, the Sobol $S_q^1$ equals the
normalized squared local sensitivity
$S_q^2\,\mathbb{V}[\theta_q]/\mathbb{V}[Y]$ and $S_q^T=S_q^1$ (no interactions).
Discrepancies between the local and global rankings are therefore *diagnostics
of nonlinearity/interaction*, not contradictions — a useful cross-check when
both are available.

---

## 17.6 Smoothness of the parameter-to-output map

Sensitivities are only meaningful if the map $\boldsymbol{\theta}\mapsto\mathbf{y}$
is (at least) continuously differentiable, and a central difference converges at
its $h^2$ rate only where $\mathbf{y}$ is smooth on $[\theta-h,\theta+h]$. Several
modeling choices in Pyrolysis.jl make the semi-discrete residual smooth — they are
precisely the smooth replacements for the classical non-smooth pyrolysis
closures, and the same smoothness keeps the structured Jacobian (§15) exact. Each
is verified against the source below. They make the *residual* smooth; a
configuration that `solve` warns about can still put a kink in it (§17.7.7).

**(1) Smooth temperature gates.** The hard Heaviside window
$H(T-T_{\min})H(T_{\max}-T)$ would make $\partial r/\partial T$ a delta function
at the gate edges, destroying the integrator's Newton Jacobian and making an
output jump wherever a perturbation moves a cell across a gate edge. The code
replaces it with the product of tanh ramps

$$
r \;\mathrel{*}=\; \tfrac12\big(1+\tanh((T-T_{\min})/w_{\text{gate}})\big)\,
                   \tfrac12\big(1+\tanh((T_{\max}-T)/w_{\text{gate}})\big),
$$

applied only when the bounds are finite/positive (`src/Physics/kinetics.jl`), with
the per-reaction gate width $w_{\text{gate}}$ (default 1 K, §6.2.2).
This is $C^\infty$ in $T$, so all temperature sensitivities of the rate are
finite. The cost is a leakage of reactivity of about $10\,w_{\text{gate}}$ outside
the nominal window (the ramp is not unit-step), discussed in §6.2; for sensitivity
purposes the benefit — globally finite $\partial r/\partial T$ — outweighs it.

**(2) Depletion limiter.** Near zero reactant the bare Arrhenius rate is stiff
and its parameter sensitivities blow up. Each reactant multiplies the rate by
the smooth factor

$$
\tanh(\xi/\xi_{\text{thr}}),
$$

with default threshold $\xi_{\text{thr}}=1.0\ \text{kg·m}^{-3}$
(`src/Physics/kinetics.jl`; reactant clamped via $\max(\xi,0)$ first).
The rate vanishes at depletion with a bounded $\partial r/\partial\xi$ for
every reaction order $n \ge 0$, keeping the reaction-rate Jacobian — and hence
the parameter sensitivities of MLR near burnout — bounded and well-defined.
The factor is $C^\infty$ in $\xi > 0$ with no threshold branch; the only
non-smooth point left is the $\max(\xi,0)$ clip exactly at $\xi=0$ (§17.7).

**(3) Smooth flux limiters and lateral-shrinkage derivative.** The TVD limiters
(Van Leer is the ALE default) are continuous, and the $\chi$-driven radial
strain rate uses an *automatic-differentiation* derivative of the
lateral-shrinkage law,
$\theta_A=\dfrac{\mathrm{law}'(\bar\chi)}{\mathrm{law}(\bar\chi)}\dfrac{\mathrm{d}\bar\chi}{\mathrm{d}t}$
with `dlaw_dχ̄ = ForwardDiff.derivative(law, χ̄)` and a guard
`A_ratio > eps()` (`src/Residual/residual.jl`). Because $\mathrm{law}'$
is taken by AD, *any differentiable* law gives an exact strain rate and
Jacobian — but a non-smooth law breaks this (see §17.7).

**(4) Harmonic-mean face averaging and tanh-based smoothing elsewhere.** The
distance-weighted harmonic mean used for face $k_f,\lambda_f,\kappa_f$ (§5.12)
and the smooth (rather than `min`/`max`) constructions in the gates keep the
discrete operators differentiable. The few genuinely non-smooth operations that
remain — `clamp` on porosity/temperature, `max` in the depletion clamp, and the
state rewrite at park events — are localized and are exactly the pitfalls
catalogued next.

> **Implementation note.** The smoothness-critical kernels are
> `reaction_rate` (`src/Physics/kinetics.jl`), the limiter library
> (`src/Discretization/flux_limiters.jl`), and `_compute_theta_A`
> (`src/Residual/residual.jl`). Each was written with `tanh`/AD
> derivatives *in place of* the textbook step/cutoff/finite-difference form so
> that the residual is differentiable in the state and the parameters.

---

## 17.7 Pitfalls

Even with a smooth residual, sensitivity computations can return wrong or
meaningless numbers if certain solver features are active or the parameter point
sits on a non-smooth manifold. The pitfalls below are specific to Pyrolysis.jl.

### 17.7.1 Warm-start artifacts

`solve(...; warm_start=prev_solution)` reuses the prior solve's *final* state as
the new initial condition (`_apply_warm_start!`, `src/Solver/solve.jl`). Inside
a sensitivity loop this introduces a **hidden parameter dependence**: the
warm-start state was produced at a *different* $\boldsymbol{\theta}$, so both
perturbed solves start from the same state and the difference sees
$\partial\mathbf{u}_0/\partial\theta_q = 0$, even though the true initial
condition (had the perturbed problem been integrated from the real start) would
depend on $\boldsymbol{\theta}$. A finite difference of warm-started solves is
therefore a derivative of the *wrong map*. **Rule:** do not warm-start the solves
of a finite-difference evaluation or a sensitivity sweep; integrate each sample
from its true initial condition.

### 17.7.2 Non-smooth lateral-shrinkage laws

The radial strain rate uses $\mathrm{law}'(\bar\chi)$ via AD (§17.6, item 3). Identity
(`NoChi`) and the smooth linear laws $1\mp k\bar\chi$ are differentiable, but a
piecewise or table-defined law with kinks makes $\mathrm{law}'$ discontinuous at
the breakpoints. A sensitivity evaluated *at* a breakpoint is then ill-defined
(left/right derivatives differ), and a central difference averages the two
one-sided slopes. The guard `A_ratio > eps()` further
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
trajectory, so $\mathbf{u}(t)$ is only piecewise smooth in the parameters.
Finite-difference sensitivities are valid on horizons without park events; a
pair of perturbed solves that park at different times measures the jump, not the
derivative.

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

### 17.7.5 Solver noise in finite differences

A finite difference measures the integration error along with the derivative
(§17.2.2). On a well-posed problem that error changes smoothly with
$\boldsymbol{\theta}$ for perturbations that leave the step sequence unchanged,
and the estimate is stable over several decades of step; a step too large incurs
truncation error, largest for pointwise outputs. Confirm a step by repeating
with $10h$ and $h/10$. Disagreement means the output is not smooth in
$\boldsymbol{\theta}$ there, most often because of one of the discrete features
of §17.7.1–§17.7.3 or the degenerate gas transport of §17.7.7.

### 17.7.6 Non-differentiable guards in the kinetics and clamps

The reactant clamp $\max(\xi,0)$, the temperature clamp to $[200,3000]\ \text{K}$
in the boundary Newton solve, and the porosity `clamp(\cdot,0,1)` are all
non-smooth at their boundaries. In normal operation states stay in the interior
and these are inactive, but a sensitivity evaluated exactly *at* a clamp boundary
(e.g. a fully depleted reactant, $\xi=0$) inherits a one-sided or zero
derivative. This is rarely an issue for QoIs integrated over time but can surface
in pointwise late-time MLR sensitivities near burnout.

### 17.7.7 Degenerate gas transport

Condensed components default to a gas-transfer coefficient $\lambda=0$. Under
Fickian gas transport the cell coefficient is then $\sum_{j\in\text{gas}} v_j\lambda_j$,
proportional to the cell's gas volume fraction (§9.2), and the face coefficient — a harmonic
mean, set to zero when either side is below $10^{-20}$ — vanishes on every face
into a cell that holds no gas. The gas field ends in a sharp front. The cell
ahead of the front gains gas only once its concentration is positive, and its
concentration sits at the integration noise, below `abstol`, crossing zero.
Once seeded, its gas grows at a rate of order $10^{2}\ \text{s}^{-1}$ on the
User Guide example, so the right-hand side has a kink of large slope at that
crossing. The implicit Newton fails repeatedly on it (the integrator's steps fall
to a few milliseconds), and the time at which each cell fills — and with it the
whole trajectory — depends on the noise. Finite differences of such a trajectory
measure the noise, not a derivative (§17.2.2). `solve` warns in this
configuration (`warn_degenerate_gas_transport`); a small $\lambda$ on the
condensed components (e.g. $10^{-6}\ \text{m}^2\,\text{s}^{-1}$) removes the
front.

---

## 17.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

None of the three reference codes differentiates its solver with respect to the
model parameters, and neither does Pyrolysis.jl: all four obtain parameter
sensitivities *externally*, through driver scripts and repeated forward runs.

**ThermaKin2Ds.** ThermaKin couples a forward pyrolysis solver to inverse
analysis via external optimization (e.g. Stochastic Hill-Climbing / Shuffled
Complex Evolution over the kinetic and property parameters), evaluating the
forward model as a black box and forming sensitivities by repeated forward runs
(finite differences). Mapping notation: ThermaKin's stoichiometric coefficients
$\theta_i^j$ are our yields $\nu_{r,j}$, and its heats of reaction carry the
*opposite* sign ($h>0$ exothermic) which our $Q_{\text{rxn}}=-\sum_r h_r r_r$
converts to (overload H1, §2). The local derivative of §17.2 is the same
black-box construction, and its step-size analysis (§17.2.2) applies to either
code: an adaptive stiff integrator leaves an error in each output that the step
must dominate.

**Gpyro.** Gpyro's calibration (Lautenberger & Fernandez-Pello) is built on a
*genetic algorithm* (GA) that minimizes the mismatch between simulated and
measured MLR/temperature over a population of parameter sets — derivative-free by
design. Gpyro maps the pre-exponential $Z$ → our $A_r$ and the conversion-based
nth-order kinetics to our mass-concentration form (§6.9). A GA explores a global
landscape without gradients (robust to the $(A,E)$ degeneracy, to non-smooth
responses and to solver noise) but needs many forward solves and gives no local
sensitivity matrix. Pyrolysis.jl's shipped calibration is likewise
derivative-free (Nelder–Mead); finite-difference gradients can drive local
refinement (Gauss–Newton/L-BFGS) once a GA or sweep has found a basin, and the
variance-based Sobol harness (§17.5) supplies the global parameter ranking that
a GA does not.

**FDS (solid phase).** The FDS solid-phase pyrolysis model recovers kinetic
parameters from a reference TGA temperature and peak rate (the
$A,E$-from-$T_{\text{peak}}$ recipe, §6.9) rather than offering a built-in
sensitivity/adjoint capability; FDS's broader uncertainty treatment is the
verification-and-validation accounting of its guides, not parametric
differentiation of the solid solver.

In short, Pyrolysis.jl's sensitivity tooling is a set of recipes on top of
`solve` — central differences for local derivatives, the sweep harness and a PCE
surrogate for global indices — sharing one immutable-problem / mutable-workspace
harness. Its AD is confined to the state Jacobian (§15).

---

## 17.9 Limitations and open issues

1. **No exact parameter derivatives (open).** Local sensitivities are central
   differences of whole solves (§17.2). Exact forward or adjoint derivatives
   need the three changes of §17.3.1 — parameter-generic material types, a
   parameter-generic workspace and integrator path, and $\boldsymbol{\theta}$
   lifted into the ODE parameter vector — and then a native integration of
   (17.1) or (17.2)–(17.3) on the structured Jacobian of §15.

2. **Sensitivities across park events (open).** A park is a discrete rewrite of
   the state, so trajectories are only piecewise smooth across events and
   finite differences are valid only on horizons without one
   (`handle_depletion=true` runs that park). No sensitivity jump condition for
   the event is formulated.

3. **Warm-start in sensitivity loops (open).** Whether — and how — reusing a
   prior solution's final state as a fresh initial condition perturbs sensitivity
   estimates is unquantified. The current guidance is to avoid
   warm-starting inside sensitivity evaluations (§17.7.1).

4. **Non-smoothness residue.** A handful of guards (`max`/`clamp`, the
   `A_ratio > eps()` switch, hard temperature/
   porosity clamps) remain non-differentiable at their boundaries (§17.7.2,
   §17.7.6). They are inactive in the interior of normal operation but can
   contaminate sensitivities evaluated exactly on the boundary; a fully smooth
   reformulation of these guards is not yet in place.

5. **Parameter handles are unregistered.** Because injection is by user closure
   (`p_inject`/`modify`), there is no enforced catalogue of "supported"
   parameters or automatic check that the injected field actually reaches the
   residual; responsibility for a clean, non-mutating, $C^1$ injection rests with
   the user (§17.4.1).

6. **Global-sensitivity tooling lives in examples, not the core.** Sobol/PCE is
   implemented in `examples/sensitivity_analysis/`, not as a packaged function;
   the design and degree-selection choices (PCE order, BCa bootstrap,
   group-correlation handling) are reference recipes rather than a stable API.

Cross-references: the structured Jacobian that an exact forward or adjoint
sensitivity would reuse is §15 (Temporal Integration & Structured Jacobian); the
smooth kinetics gates and depletion limiter are §6 (Reaction Kinetics); the
mixture-property rules are §5 (Effective Properties); the lateral-shrinkage law
and $\chi$ dynamics are §10 (Volume Change & Lateral Shrinkage); the ALE geometry
is §11 (The ALE Framework); and the conservation diagnostics used to
sanity-check sensitivity-driven calibrations are §16 (Conservation
Diagnostics). Validation of sensitivity-driven fits against experiment and
reference codes is §18 (Verification & Validation).
