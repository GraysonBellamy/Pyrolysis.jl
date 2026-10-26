# 17. Sensitivity Analysis (Theory)

This chapter develops the mathematics of parametric sensitivity for the
Pyrolysis.jl solver: how a derivative of a simulation output with respect to
model parameters is obtained, how accurate it is, and how it supports inverse
parameter estimation, uncertainty quantification, and optimization. We treat
both *local* sensitivity — central finite differences of the entire `solve`
chain — and *global*, variance-based sensitivity (Sobol indices computed through
the parameter-sweep harness). Pyrolysis.jl gives sensitivities by finite
differences of whole solves; no derivative is taken through the solve. The
solve path is `Float64` throughout: the material constructors convert every
parameter to `Float64`, and the integrator's right-hand side, Jacobian and
solution containers are `Float64`-typed. Neither a `ForwardDiff.Dual` nor a
reverse-mode pullback can therefore pass from a parameter to an output
(§17.3 explains what exact derivatives would require). Neither flavor has a
dedicated API; both are short user scripts around `solve`, and the mathematics
below is what those scripts compute.

Conventions used in this chapter (binding, from §2, Nomenclature):

- $\boldsymbol{\theta}\in\mathbb{R}^{p}$ is the vector of **parameters** with
  respect to which sensitivities are taken (e.g. $A_r$, $E_r$, $h_r$, a density
  $\rho_j$, a conductivity $k_j$, or a boundary-condition coefficient). $p$ is
  the parameter count.
- $\mathbf{u}(t)\in\mathbb{R}^{N_{\text{state}}}$ is the state of the
  mass-matrix system $M\,\dot{\mathbf u}=\mathbf f(\mathbf u,t)$ the integrator
  advances (§15): per cell the component masses and the energy state (the
  temperature, or the sensible-enthalpy inventory); on a moving mesh the node
  positions and node velocities, the velocities algebraic rows; under
  `BEER_LAMBERT` the face intensities, algebraic rows; under a lateral
  shrinkage law the conversion-history pair of every cell and the
  cross-section rows. $M$ is diagonal, with zeros on the algebraic rows.
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
release rate, time to ignition, residual char fraction, burnout time. The
**sensitivity** of an output to a parameter,

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
| Local | central differences of `solve` (§17.2) | $2p+1$ solves per step | small $p$ (a handful of parameters), any number of outputs |
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
back inside $p_{\text{inject}}$. Equivalently, a parameter perturbed relatively,
$p=p_0(1+\varepsilon)$, gives the central difference
$D(h)=\big(J(h)-J(-h)\big)/(2h)$ of an objective $J(\varepsilon)$, an estimate of
$dJ/d\ln p$.

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
step sequence and the step rejections change discretely as $\theta$ moves. A
perturbation small enough to change none of those decisions leaves $\epsilon$
smooth; one that changes a decision makes the noise term of order
$\lvert\epsilon\rvert/h$. Balancing that bound against truncation gives the
classical optimum $h^\ast\sim\big(3\lvert\epsilon\rvert/\lvert y'''\rvert\big)^{1/3}$.

On the charring-slab example of User Guide §13.2 (40 cells, 50 kW m⁻²,
`reltol = 1e-7`, `abstol = 1e-9`, the default integrator) the estimates are
stable over several decades of step: the derivatives of the peak and 150-s
surface temperatures, the mean surface temperature and the final mass with
respect to $(\log_{10}A, E, h)$ at $\eta=10^{-5}$–$10^{-7}$ agree with those at
$\eta=10^{-4}$ to $3\times10^{-4}$ relative or better (the final mass to
$10^{-5}$). Truncation error appears at larger steps — up to 0.4 % at
$\eta=10^{-3}$ and 28 % at $\eta=10^{-2}$, largest for the 150-s temperature,
which moves with the reaction front.

The same example with the condensed components at their default $\lambda=0$
behaves very differently (§17.7.7): the solve accepts 2,016 steps and rejects
379 (295 and 23 with $\lambda=10^{-6}$), a $10^{-6}$ relative change in $A$ moves
$T_s$ by up to 0.29 K ($8\times10^{-6}$ K with $\lambda=10^{-6}$), and the
derivative of the peak surface temperature with respect to $\log_{10}A$ reads
$-1.08$, $-0.83$, $-3.1$ and $+48.8$ K at $\eta=10^{-3}$, $10^{-4}$, $10^{-5}$ and
$10^{-6}$, against $-1.123$ K at every one of those steps with $\lambda=10^{-6}$.

**A convergence check.** A quick confirmation of a step $\eta$ is to repeat with
$10\eta$ and $\eta/10$ and treat disagreement as a sign that $y$ is not smooth in
$\boldsymbol\theta$ at this point (§17.7). A derivative that has to be relied on
is qualified by a convergence check on the step and the tolerance together.
Perturb each parameter relatively by $\pm h$ at $h = 1, 2, 4$ and $8\times10^{-3}$
and at two tolerances (`reltol = 1e-7` and `1e-8`, `abstol` scaled with
`reltol`), form the central differences $D(h)$ and the Richardson extrapolants
of successive steps,

$$
R(h) = \frac{4\,D(h)-D(2h)}{3},
$$

which remove the $h^2$ truncation term, and take the value where the
extrapolants of two successive pairs agree within 5 % and the two tolerances
agree at those steps. A step too large carries the objective's curvature; a
step too small carries the solve's error, or the objective's own roughness
where it samples an event, such as a deactivation, at the data's times. For
the V-BO-01 burnout deck (§18; 60 PMMA cells, the mean square misfit of the
mass-loss rate against Gpyro V2) the check gives $dJ/d\ln\alpha=-24.65$,
$dJ/d\ln q=957.7$ and $dJ/d\ln A=32.57$, with steps of 1 to $4\times10^{-3}$,
2 to $8\times10^{-3}$ and 1 to $2\times10^{-3}$ usable within 5 % at `reltol`
$10^{-6}$ to $10^{-8}$; a step of $2\times10^{-3}$ at `reltol = 1e-7` serves all
three. These numbers are that case's; the check is repeated for another deck
or objective.

**Steps through an event.** Where the solution passes an event, a located
event's time is a smooth function of the parameters (it is the first root of
the threshold margin, §14), so the derivative exists. An objective that
samples the solution at fixed times sees the event's burst move across its
sample times, which adds a roughness of the period of one sample interval in
the event time. On V-BO-01 the deactivation moves 0.46 s per $10^{-3}$ of flux
against the reference's 1 s samples, and the $10^{-3}$ difference departs by
6 % at every tolerance, the misfit before the burnout staying smooth to
$10^{-5}$; the convergence check's usable steps exclude such a step. The
derivative through the solve stays outside the package (§17.3).

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
the **forward sensitivity equation** obtained by differentiating the
mass-matrix system $M\dot{\mathbf{u}}=\mathbf{f}(\mathbf{u},\boldsymbol{\theta},t)$:

$$
M\,\dot{\mathbf{s}}_q
=\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\,\mathbf{s}_q
+\frac{\partial\mathbf{f}}{\partial\theta_q},
\qquad
\mathbf{s}_q(0)=\frac{\partial\mathbf{u}_0}{\partial\theta_q},
\tag{17.1}
$$

whose algebraic rows (zero rows of $M$) are constraints the sensitivity must
satisfy at every time, its initial value included. Here
$\partial\mathbf{f}/\partial\mathbf{u}$ is precisely the colored Jacobian of
§15. Pyrolysis.jl does **not** integrate (17.1); the central difference
approximates $\mathbf{s}_q$ by the difference of two primal trajectories,
including the contribution of $\partial\mathbf{u}_0/\partial\theta_q$ when a
parameter appears in the initial condition (an initial density, which also sets
the cells' reference condensed masses and with them the mass tolerances and the
deactivation threshold), because $p_{\text{inject}}$ rebuilds the problem and
`setup` rebuilds the initial state from it. Equation (17.1) is the object a
native implementation would integrate (§17.3.2).

### 17.2.4 Cost and parallelism

A central-difference Jacobian costs $2p$ perturbed solves plus the base solve,
independent of the number of outputs $m$: the full $m\times p$ matrix $S$ comes
from the same solves. That suits Gauss–Newton, whose normal equations need $S$
itself, and a handful of parameters. The convergence check of §17.2.2 costs
$8p+1$ solves per tolerance. The perturbed solves are independent;
`parameter_sweep` (§17.5.3) runs a list of perturbed parameter vectors with
$p_{\text{inject}}$ as its `modify` and distributes them over `Distributed`
workers (User Guide §13.5). A gradient of a scalar loss costs the same $2p$
solves, so for large $p$ a derivative-free or global method is usually the
better investment.

---

## 17.3 Why the solve is not differentiated automatically

### 17.3.1 What stops dual numbers and adjoints

The residual kernels are scalar-type-generic — they run on `ForwardDiff.Dual`
states for the colored Jacobian (§15) — but the path from a *parameter* to an
output is not:

1. **Construction.** Scalar properties go through
   `to_property(x::Real) = ConstantProperty(Float64(x))`
   (`src/Materials/property_functions.jl`), and `Reaction` stores $A$ and $E$ as
   `Float64` (`src/Materials/reactions.jl`). A dual-valued parameter fails with
   `MethodError: no method matching Float64(::ForwardDiff.Dual{…})` before any
   solve starts.
2. **Set-up problem and integrator.** `setup` builds the solve's definition from
   those `Float64` parameters; its workspace holds `Float64` buffers and
   dual-number twins for the state's derivative only. The right-hand side and
   the Jacobian reach the integrator behind fixed-signature wrappers for
   `Vector{Float64}` states and `Float64` matrices (`ErasedRHS`,
   `ErasedJacobian`), so that the integrator compiles once for every problem
   (§15); `PyrolysisSolution` stores its traces as `Float64` arrays.
3. **No parameter vector.** The `ODEProblem` that `integrate` builds carries no
   parameter vector (`p = NullParameters()`); $\boldsymbol{\theta}$ enters only
   through the set-up problem, where an adjoint sensitivity algorithm cannot see
   it.

### 17.3.2 What exact derivatives would require

Exact forward derivatives would integrate (17.1) alongside the primal, reusing
the colored Jacobian $\partial\mathbf{f}/\partial\mathbf{u}$ of §15 and
evaluating $\partial\mathbf{f}/\partial\theta_q$ by forward-mode AD of the
residual with respect to $\boldsymbol{\theta}$ alone. For a scalar loss
$\ell=\int_0^{t_{\text{end}}}g(\mathbf{u},\boldsymbol{\theta},t)\,\mathrm{d}t+\Phi(\mathbf{u}(t_{\text{end}}),\boldsymbol{\theta})$
and many parameters, the continuous adjoint is cheaper: the costate
$\boldsymbol{\lambda}(t)$ solves the backward equation

$$
M^{\mathsf T}\dot{\boldsymbol{\lambda}}
=-\left(\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\right)^{\!\mathsf T}\boldsymbol{\lambda}
-\left(\frac{\partial g}{\partial\mathbf{u}}\right)^{\!\mathsf T},
\tag{17.2}
$$

with terminal values $\boldsymbol\lambda(t_{\text{end}})$ set from
$\partial\Phi/\partial\mathbf u$ on the differential rows and consistent with the
algebraic rows, after which the gradient is the quadrature

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
transpose of the colored Jacobian as its operator. Either route needs the
parameters lifted into a parameter vector the residual reads (item 3 above) and
a residual that is generic in the parameter element type (item 1). Both also
meet the deactivation events (§17.7.3): an event is a discrete rewrite of the
state at a located time — inventories transferred to the receiver, rows
masked — and carrying a sensitivity through it needs a jump condition (the
event time's sensitivity from its root condition and the derivative of the
transfer map), which is not formulated.

---

## 17.4 The `p_inject` wrapper pattern and parameter handles

### 17.4.1 Why wrap rather than rewrite

The residual reads its parameters from the definition `setup` builds from the
problem — the material, the boundary conditions, the mesh — and has no
parameter argument. Rather than thread a parameter vector
$\boldsymbol{\theta}$ through the definition and every kernel, sensitivities are
obtained by re-running *problem construction* followed by the solve at each
parameter point. The user supplies a closure

$$
p_{\text{inject}}:\ (\text{base\_problem},\boldsymbol{\theta})\ \longmapsto\ \text{problem}',
$$

that materializes a *fresh* `PyrolysisProblem` with $\boldsymbol{\theta}$ baked
into the relevant fields, and the finite difference (or the Sobol sampler)
evaluates $\boldsymbol{\theta}\mapsto\texttt{solve}(p_{\text{inject}}(\dots))$.

This has three consequences worth stating precisely:

1. **No invasive refactor.** The residual, the kernels, the colored Jacobian
   and the state layout are untouched. Problem construction (building a
   `Material` and `BoundaryConditionSet`) and `setup` are re-run at every
   evaluation; both are cheap next to the solve.
2. **Any constructible parameter can be perturbed.** Because $p_{\text{inject}}$
   builds the problem, *anything* that flows into a `Material`,
   `BoundaryConditionSet`, mesh, or initial condition is a valid handle — there
   is no fixed parameter registry.
3. **The base problem must not be mutated.** `p_inject` returns a *new*
   problem; mutating `base_problem` in place would make one evaluation leak into
   the next (the same contract the parameter-sweep harness enforces, §17.5).
   The new problem may share the base problem's mesh, which a solve never
   changes.

### 17.4.2 Parameter handles

Common $\boldsymbol{\theta}$ entries and how they enter the model:

- **Kinetic parameters** $A_r$, $E_r$. Rebuild the reaction(s) with new
  pre-exponential/activation energy and assemble a new `Material`. The rate
  $r=A\,e^{-E/(R_gT)}\prod_j\xi_j^{n_{r,j}}$ (times the depletion limiter's
  factor and the temperature gates, §6) is analytic in both $A$ (linear) and
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
  (and hence $\chi$, §10) and the cells' initial condensed masses, which scale
  the mass rows' absolute tolerances and set the deactivation threshold
  (§14), so density is a particularly broad-reaching handle.
- **Boundary-condition parameters** — incident flux level, convective
  coefficient $h_{\text{conv}}$, ambient $T_\infty$, view factor $F$,
  emissivity in $q_{\text{rad}}=F\varepsilon\sigma(T_\infty^4-T_s^4)$. These
  enter the surface energy balance solved for $T_s$ (§12.2).

### 17.4.3 Worked structure of a local-sensitivity evaluation

```julia
using Pyrolysis, FiniteDiff

# The base problem: a charring slab under 50 kW/m² (User Guide §13.2).
virgin = SolidComponent(:virgin; ρ = 1190.0, φ=.3, c = 1500.0, k = 0.20, λ = 1e-6)
char = SolidComponent(:char; ρ = 200.0, φ=.3, c = 1200.0, k = 0.10, λ = 1e-6)
gas = GaseousComponent(:gas; M = 0.100, c = 1500.0, k = 0.02,)
rxn0 = Reaction(:decomp, 1 => (2, 3); A = 1e13, E = 1.88e5, h = 8.7e5, yields = (0.2, 0.8))
base_problem = PyrolysisProblem(;
    mesh = create_uniform_mesh(10.0e-3, 40),
    material = Material(name = :charring_polymer, components = (virgin, char, gas),
                        reactions = (rxn0,)),
    bc_set = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0)
                         + ConvectiveBC(h = 10.0, T_inf = 300.0),
                   bottom = AdiabaticBC()),
        mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
    ),
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0, z -> 0.0],
    tspan = (0.0, 300.0),
    output_times = collect(range(0.0, 300.0, length = 151)),
)

# θ = [log10(A), E, h]: a fresh problem with θ in its reaction.
function p_inject(prob, θ)
    m = prob.material
    rxn = Reaction(:decomp, 1 => (2, 3); A = 10.0^θ[1], E = θ[2], h = θ[3],
                   yields = (0.2, 0.8))
    mat = Material(; name = m.name, components = m.components, reactions = (rxn,),
                   mixing = m.mixing, lateral_shrinkage_law = m.lateral_shrinkage_law,
                   depletion_limiter = m.depletion_limiter, gas_viscosity = m.gas_viscosity, gas_transport = m.gas_transport)
    return PyrolysisProblem(; mesh = prob.mesh, material = mat, bc_set = prob.bc_set,
                            T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
                            tspan = prob.tspan, dt_initial = prob.dt_initial,
                            output_times = prob.output_times)
end

θ0 = [13.0, 1.88e5, 8.7e5]
y(θ) = [solve(p_inject(base_problem, θ); abstol = 1e-9, reltol = 1e-7).extras.total_mass[end]]

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
(User Guide Chapter 12). Its contract dovetails with the sampling above:

$$
\texttt{modify}(\text{base\_problem},\,\mathbf{p}_s)\ \longmapsto\ \text{problem}_s,
\qquad
\texttt{extract}(\text{solution}_s)\ \longmapsto\ \mathbf{y}_s,
$$

for each sample point $\mathbf{p}_s$ in the design `grid`. `modify` plays the
same role as `p_inject` in §17.4 — it must build a *fresh* problem and must not
mutate the base. The harness dispatches the (embarrassingly parallel) solves
over `Distributed.pmap` when workers are available; each solve sets up its own
problem, workspace and Jacobian, so nothing mutable is shared between points:

```julia
work = point -> extract(solve(modify(problem, point); solve_kwargs...))
points = collect(grid)
results = parallel ? pmap(work, points) : map(work, points)   # parallel = nworkers() > 1
```

(`src/Integration/parameter_sweep.jl`.) The `extract` closure should return small
NamedTuples of scalar QoIs (peak $T_s$, integrated MLR, time-to-peak) to
minimize worker→head transfer; those scalars are exactly the $Y$ samples the PCE
surrogate consumes. Thus the *same* solve and the *same* rebuild-do-not-mutate
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
the smooth replacements for the classical non-smooth pyrolysis closures, and the
same smoothness gives the colored Jacobian (§15) a finite, continuous value at
every state. Each is verified against the source below. They make the
*residual* smooth; a configuration that `solve` warns about can still put a
kink in the trajectory (§17.7.7).

**(1) Smooth temperature gates.** The hard Heaviside window
$H(T-T_{\min})H(T_{\max}-T)$ would make $\partial r/\partial T$ a delta function
at the gate edges, destroying the integrator's Jacobian and making an output
jump wherever a perturbation moves a cell across a gate edge. The code
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
the one-sided factor $g(\xi/\xi_{\text{thr}})$ of the `DepletionLimiter`
(`src/Materials/reactions.jl`, `one_sided_factor` in `src/Physics/kinetics.jl`;
default threshold $\xi_{\text{thr}}=1.0\ \text{kg·m}^{-3}$): exactly zero at and
below zero concentration, a cubic onset of zero value and slope at zero, joined
with matching value and slope onto $\tanh(\xi/\xi_\text{thr})$ at the join
(default one threshold), and $\tanh$ above it. The reaction-order power is taken
of the finite base $\xi_b\,\mathrm{softplus}(\xi/\xi_b)$ (scale
$\xi_b=10^{-8}$ kg·m⁻³), which equals $\xi$ to within $\xi_b$ and is finite for
every real $\xi$. The rate is therefore C¹ through zero concentration for every
reaction order $n\ge0$, with a bounded $\partial r/\partial\xi$, keeping the
reaction-rate Jacobian — and hence the parameter sensitivities of MLR near
burnout — bounded and well-defined; far above the threshold the bare rate law is
recovered.

**(3) Limited face values and the lateral-shrinkage rate.** The face values of
the relative flux on a moving mesh are limited by a smooth (van Albada)
limiter (`src/Kernels/relative_flux.jl`). Under a lateral shrinkage law the
area follows $A=A_0\,\mathrm{law}(\bar\chi)$, and the lateral strain rate is the
automatic-differentiation derivative of $\log\mathrm{law}(\bar\chi)$ along the
column's conversion rates (`ForwardDiff.derivative` in
`src/Closures/lateral.jl`), with $\bar\chi$ smoothly clamped into
$[-10^{-3},\,1+10^{-3}]$ before the law. The law is checked positive and finite on
that interval when the problem is set up, so the area never vanishes and no
guard switches the rate off. Because $\mathrm{law}'$ is taken by AD, *any
differentiable* law gives an exact strain rate and Jacobian — but a non-smooth
law breaks this (see §17.7.2).

**(4) Physical branches and trial guards.** Conductivity uses a guarded
resistance mean. Gas transport uses total-flux interface roots, explicit
vacuum branches and deterministic donor selection. Its implicit derivative
requires a nonsingular qualified root; support changes and donor switches
need separate sensitivity checks (§17.7.7). Trial guards on temperature,
masses, volumes and capacities are C¹ extensions, with accepted-state checks
ensuring the physical solution remains in their identity region. Gas guard
scales track current pore volume and temperature. Integrator decisions and
deactivation events remain discrete operations.

> **Implementation note.** The smoothness-critical kernels are
> `reaction_rate` and `one_sided_factor` (`src/Physics/kinetics.jl`), the
> guards (`src/Materials/guards.jl`), `guarded_harmonic_mean`
> (`src/Properties/effective_properties.jl`), the limited face values of
> `src/Kernels/relative_flux.jl`, and the lateral rows (`src/Closures/lateral.jl`).
> Each is written with `tanh`, polynomial joins or AD derivatives *in place of*
> the textbook step, cutoff or finite-difference form, so that the residual is
> differentiable in the state and the parameters.

---

## 17.7 Pitfalls

Even with a smooth residual, sensitivity computations can return wrong or
meaningless numbers if certain solver features are active or the parameter point
sits on a non-smooth manifold. The pitfalls below are specific to Pyrolysis.jl.

### 17.7.1 Solve options that differ between the perturbed solves

A finite difference of `solve` is the derivative of the *discretized* model:
the mesh, the tolerances and every option that changes the solution are part
of the map. Every solve of one difference — and of one convergence check at a
given tolerance — must be made with identical options: `energy_state`,
`abstol`, `reltol`, `use_ale`, `radiation_model`, `deactivation`, `threshold`,
`band`, `located_events`, `integrator_policy`, `linear_solver`, `dtmax`,
`saveat`, and the problem's `output_times`. An option that differs between the
$+h$ and the $-h$ solve adds its own effect divided by $2h$ to the estimate.
**Rule:** pass the options through one named tuple to every solve (User Guide
§13.3), and keep `tspan`, `output_times` and the mesh literals of the base
problem rather than functions of $\boldsymbol\theta$.

### 17.7.2 Non-smooth lateral-shrinkage laws

The lateral strain rate uses $\mathrm{law}'(\bar\chi)$ via AD (§17.6, item 3).
Identity and the smooth linear laws $1\mp k\bar\chi$ are differentiable, but a
piecewise or table-defined law with kinks makes $\mathrm{law}'$ discontinuous at
the breakpoints. A sensitivity evaluated *at* a breakpoint is then ill-defined
(left/right derivatives differ), and a central difference averages the two
one-sided slopes. **Rule:** use $C^1$ shrinkage laws for sensitivity studies,
and avoid evaluating at $\bar\chi$ values where the law has a corner.

### 17.7.3 Deactivation events

Burnout is handled by deactivation events (§14): at the first time a live
cell's condensed mass reaches its threshold, located inside the step, the cell
— on a moving mesh with the cells of the near-threshold band outward of it —
is merged into the cell beneath and its rows are masked. The located event's
time is a smooth function of the parameters, so the trajectory after the event
moves smoothly with $\boldsymbol\theta$ as long as the event takes the same
cells. Two things break that smoothness:

1. **Sampled objectives.** An objective that samples the solution at fixed
   times sees the event's burst move across its sample times (§17.2.2); its
   finite differences carry a roughness of the period of one sample interval
   in the event time, which the convergence check's usable steps exclude.
2. **A change of the event's decision.** Which cells one invocation takes is
   decided by threshold tests on the state; a perturbation that moves a cell
   into or out of the band, or splits one invocation into two, changes the
   merged state by a finite amount. A pair of perturbed solves on the two sides
   of such a change measures the jump, not the derivative. The number of
   invocations is in `sol.statistics.invocations` and every event in
   `sol.events`; comparing them between the perturbed solves detects the case.

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
with $10h$ and $h/10$, or qualify it with the convergence check of §17.2.2.
Disagreement means the output is not smooth in $\boldsymbol{\theta}$ there, most
often because of one of the discrete features of §17.7.1–§17.7.3 or the
degenerate gas transport of §17.7.7.

### 17.7.6 Upwind switches and the C¹ guards

The Darcy advective flux takes the upwind cell's concentration on the sign of
the face velocity, and the relative flux on a moving mesh its upwind stencil on
the sign of the mesh correction (§9, §11); the flux is continuous across the
switch but its derivative is not, so a state with a face velocity at zero is a
kink of the residual. In normal operation the gas flows one way through each
face and the switch is inactive; a flow that reverses at a face, or a
correction that changes sign there, puts the kink on the trajectory. The guards
of §17.6, item 4, are C¹: their second derivative jumps at the joints, which
degrades the central difference's $h^2$ truncation rate on an interval that
straddles a joint but leaves the derivative continuous. They are the identity
on the physical states, so this is rarely an issue for QoIs integrated over
time, but it can surface in pointwise late-time sensitivities near burnout,
where a cell's masses approach their mass scale.

### 17.7.7 Transport support and implicit derivatives

Closed wall or Darcy channels are deliberate constitutive choices. Their
combination at an interface is resolved by matching total species flux.
When both paths are degenerate, composition may be unobservable and its
continuation must be stated; a derivative through a support change is not
implied by a successful residual evaluation.

General local roots require a nonsingular, numerically resolved Jacobian.
The implementation restores first state derivatives by the implicit function
theorem at that root. Vacuum and inactive channels use explicit continuations;
donor switches select a deterministic primal branch. Higher derivatives or
parameter derivatives through a singular channel transition require separate
qualification. Endpoint and capillary regularizations are constitutive inputs
whose scales must be varied in solution studies.

---

## 17.8 Comparison to ThermaKin2Ds, Gpyro, and FDS

The sources support specific calibration workflows, not a universal claim
that none of these codes can provide parameter derivatives.

**ThermaKin2Ds.** TK describes the forward solver and examples using
experimentally determined or fitted material parameters. It does not specify
a universal optimizer or sensitivity algorithm. A statement about hill
climbing, shuffled complex evolution or finite differences must identify the
particular inverse-study implementation; these are not inferred capabilities
of the forward executable.

**Gpyro.** GP's property-estimation chapter documents a genetic algorithm
using repeated forward evaluations. That procedure is derivative-free, not
a finite-difference sensitivity calculation. It does not guarantee a unique
fit or robustness to parameter non-identifiability, and does not establish
that every other Gpyro calibration workflow uses the same algorithm.

**FDS.** Its reference-temperature/rate prescription estimates Arrhenius
parameters (§6.9.3). That feature is distinct from whole-solve parameter
sensitivities. It is not evidence for a blanket absence of sensitivity or
optimization workflows around FDS.

**Pyrolysis.jl.** The supported workflow in this repository is qualified
central differences of complete solves (§17.2), parameter sweeps and
variance-based analysis. The inverse white-pine example uses Nelder–Mead.
Colored AD of the state residual supplies the integration Jacobian; it does
not provide exact whole-solve parameter gradients. These statements follow
the source paths and examples cited earlier in this chapter.

---

## 17.9 Limitations and open issues

1. **No exact parameter derivatives (open).** Local sensitivities are central
   differences of whole solves (§17.2). Exact forward or adjoint derivatives
   need the three changes of §17.3.1 — parameter-generic material types, a
   parameter-generic set-up and integrator path, and $\boldsymbol{\theta}$
   lifted into the ODE parameter vector — and then a native integration of
   (17.1) or (17.2)–(17.3) on the colored Jacobian of §15.

2. **Sensitivities across deactivation events (open).** A located event's time
   is smooth in the parameters, so derivatives exist across events that take
   the same cells; objectives sampled at fixed times add a roughness of one
   sample interval in the event time, and a change of which cells an event takes
   is a jump (§17.7.3). No sensitivity jump condition for the event is
   formulated for a native derivative.

3. **The step must be qualified per case.** The usable steps depend on the
   deck, the objective and the tolerance (§17.2.2); the convergence check is a
   recipe applied by the user, and its numbers for one case do not transfer to
   another.

4. **Non-smoothness residue.** The upwind switches of the Darcy and relative
   fluxes are kinks of the residual at zero face velocity, and the guards are
   C¹ rather than smooth (§17.7.6). They are inactive in the interior of normal
   operation but can contaminate sensitivities evaluated where a flow reverses
   or a cell approaches empty.

5. **Parameter handles are unregistered.** Because injection is by user closure
   (`p_inject`/`modify`), there is no enforced catalogue of "supported"
   parameters or automatic check that the injected field actually reaches the
   residual; responsibility for a clean, non-mutating, $C^1$ injection rests with
   the user (§17.4.1).

6. **Global-sensitivity tooling lives in examples, not the core.** Sobol/PCE is
   implemented in `examples/sensitivity_analysis/`, not as a packaged function;
   the design and degree-selection choices (PCE order, BCa bootstrap,
   group-correlation handling) are reference recipes rather than a stable API.

Cross-references: the colored Jacobian that an exact forward or adjoint
sensitivity would reuse, and the mass-matrix system it belongs to, are §15
(Temporal Integration); the smooth kinetics gates and depletion limiter are §6
(Reaction Kinetics); the mixture-property rules and face means are §5 (Effective
Properties); the lateral-shrinkage law and $\chi$ dynamics are §10 (Volume Change
and Lateral Shrinkage); the moving mesh is §11 (The ALE Framework); the
deactivation events are §14; and the conservation ledgers used to sanity-check
sensitivity-driven calibrations are §16 (Conservation Ledgers and the Audit).
Validation of sensitivity-driven fits against experiment and reference codes is
§18 (Verification and Validation).
