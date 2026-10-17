# 13. Sensitivity Analysis

*Sensitivity analysis* asks how the outputs of a pyrolysis simulation change when
its inputs change. Pyrolysis.jl supports two complementary flavors: **local**
sensitivities (`∂output/∂θ` at one parameter point, by finite differences of
whole solves) and **global / variance-based** (Sobol) sensitivities built on the
parameter-sweep harness from Chapter 12. Neither has a dedicated API; both are
short scripts around `solve`. This chapter shows how to write them, what can
legitimately be a parameter, how to wire your parameters into a fresh problem,
how to choose a finite-difference step the solver's own noise does not swamp,
and the pitfalls that make a sensitivity ill-defined. For the mathematics — the
finite-difference error balance and the variance decomposition — see Technical
Reference §17.

> **Coordinate and symbol conventions.** `z = 0` is the bottom/substrate
> (boundary id 2, tag `:bottom`); `z = L` is the top/exposed surface (boundary id
> 1, tag `:top`); heat enters from the top. The Arrhenius parameters are the
> pre-exponential `A_i` [s⁻¹] (or [m³·kg⁻¹·s⁻¹] for higher-order reactions), the
> activation energy `E_i` [J·mol⁻¹], and the heat of reaction `h` [J·kg⁻¹] with
> the **`h > 0` endothermic** storage convention (Technical Reference §2 and §6).
> Mass concentration `ξ_j` [kg·m⁻³] is the primary species state.

---

## 13.1 Overview: when to use which method

| Method | How | Returns | Cost | Best when |
|---|---|---|---|---|
| **Local** | central differences of `solve` with `FiniteDiff` (§13.3) | output + Jacobian `∂out/∂θ` | `2p + 1` solves for `p` parameters | a handful of parameters, one or many outputs, gradients for calibration |
| **Global (Sobol)** | `parameter_sweep` + a Sobol design (§13.8) | variance indices per parameter/group | `O(N_samples)` solves | uncertainty quantification, importance ranking over wide ranges |

Local sensitivities linearize about a single parameter point `θ`. Sobol is
**global**: it explores a distribution of inputs and attributes output variance
to inputs, but it gives no derivative.

Neither needs an extension. The local recipe uses `FiniteDiff`, which you add to
your own environment (`Pkg.add("FiniteDiff")`); the Sobol path is
`parameter_sweep` plus your choice of sampling package.

Automatic differentiation through `solve` is not available. The material
constructors convert every number to `Float64` and the integrator is
`Float64`-typed, so a `θ` carrying `ForwardDiff.Dual` numbers fails at
construction (`MethodError: no method matching Float64(::ForwardDiff.Dual{…})`)
and reverse-mode tools have no parameter vector to differentiate (Technical
Reference §17.3).

---

## 13.2 The `p_inject` pattern: turning a vector into a problem

Both methods share one idea: a **closure** that rebuilds a fresh
`PyrolysisProblem` from a **base problem** and a **parameter vector `θ`**:

```julia
p_inject(base_problem, θ) -> PyrolysisProblem
```

The local recipe evaluates `θ -> output(solve(p_inject(base, θ)))` at perturbed
copies of `θ`; a Sobol study evaluates it at sampled `θ`. You never touch the
residual.

> **Pyrolysis types are immutable: rebuild, do not mutate.** There is no
> `Material(existing; ...)` or `PyrolysisProblem(existing; ...)` copy
> constructor. `p_inject` must construct the components, reactions, material,
> mesh, BCs and problem **from scratch**, substituting the entries of `θ` where
> they belong — exactly like the `modify` closure for a parameter sweep
> (Chapter 12), and exactly like the `build_and_solve(θ)` function in
> `examples/sensitivity_analysis/sensitivity_douglas_fir.jl`. The base problem
> is only there to supply the parts of the configuration you are *not* varying.

A minimal, correct `p_inject` for a single charring reaction whose
pre-exponential `A`, activation energy `E`, and heat of reaction `h` we want to
perturb (with `A` passed as `log10(A)`, see the rules below):

```julia
using Pyrolysis

# θ = [log10(A), E, h]
function p_inject(base, θ)
    A, E, h = 10.0^θ[1], θ[2], θ[3]

    # λ > 0 on the condensed components keeps gas transport non-degenerate (§13.9)
    virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20, λ = 1e-6)
    char   = SolidComponent(:char;   ρ = 200.0,  c = 1200.0, k = 0.10, λ = 1e-6)
    gas    = GaseousComponent(:gas;   M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)

    rxn = Reaction(:decomp, 1 => (2, 3);
                   A = A, E = E, h = h,        # ← the perturbed parameters
                   yields = (0.2, 0.8), n = 1.0)

    material = Material(
        name       = :charring_polymer,
        components = (virgin, char, gas),
        reactions  = (rxn,),
    )

    mesh = create_uniform_mesh(10.0e-3, 40, Val(3); cross_section_area = 1.0)

    bc_set = BoundaryConditionSet(
        thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0)
                            + ConvectiveBC(h = 10.0, T_inf = 300.0),
                   bottom = AdiabaticBC()),
        mass    = (top = DiffusiveMassBC(Inf), bottom = ImpermeableBC()),
    )

    return PyrolysisProblem(;
        mesh     = mesh,
        material = material,
        bc_set   = bc_set,
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1190.0, z -> 0.0, z -> 0.0],
        tspan    = (0.0, 300.0),
        output_times = collect(range(0.0, 300.0, length = 151)),
    )
end
```

A few rules that keep `p_inject` consistent:

- **Everything not in `θ` is fixed.** Read fixed values from `base` if it is
  convenient, or just hard-code them as above.
- **Keep the time grid and the mesh fixed.** `tspan`, `output_times`, mesh
  thickness and `n_cells` should be plain literals, never functions of `θ`, so
  every perturbed solve reports the same outputs at the same times. You perturb
  physical parameters, not the time grid or the mesh count.
- **Scale the parameters.** `A ~ 10¹³`, `E ~ 10⁵` and `h ~ 10⁶` differ by many
  orders of magnitude, and a finite-difference step is relative to each entry.
  Pass `log10(A)` rather than `A` and apply the inverse transform inside
  `p_inject`, as above; the same scaling helps an optimizer (§13.7).
- **`yields` must sum to 1.0** for mass conservation (Technical Reference §6.5);
  `Reaction(...)` validates this at construction. If you make a yield a parameter
  you must rebuild the complementary yields so they still sum to 1.

### What can be a parameter `θ`

Anything that enters a constructor as a number can be a parameter. Common, useful
handles:

| Parameter | Where it lives | Notes |
|---|---|---|
| `A_i`, `E_i` | `Reaction(...; A=, E=)` | the dominant kinetic sensitivities; pass `log10(A)` |
| `h` (heat of reaction) | `Reaction(...; h=)` | scalar, tuple `(a,b)` → linear, or a closure |
| `ρ_j`, `c_p`, `k_j` | `SolidComponent(...; ρ=, c=, k=)` | pass a number, or a `CallableProperty` closing over `θ` |
| `λ_j`, `κ_j` | component `λ=`, `κ=` | gas-transfer / permeability for Darcy mode |
| `ε_j`, `α_j` | component `ε=`, `α=` | emissivity / mass-absorption |
| BC parameters | `HeatFluxBC`, `ConvectiveBC(h=,T_inf=)`, `RadiativeFluxBC(flux=)`, … | external flux, `h_conv`, `T_∞` |
| mixing `β` | `MixingSpec(WEIGHTED, β)` | conductivity / transport blend weight |

For a temperature-dependent property whose coefficients are parameters, build the
property as a closure over `θ`:

```julia
# θ = [a, b]; c_p(T) = a + b*T
cp = CallableProperty(T -> θ[1] + θ[2] * T)
virgin = SolidComponent(:virgin; ρ = 1190.0, c = cp, k = 0.20)
```

`Reaction` accepts a scalar `h` (a constant), a tuple `h = (a, b)` (linear in `T`,
`a + b*T`), or a property object such as `StateDependentProperty((T, ξ) -> ...)`
for moisture-dependent heats of sorption (Technical Reference §6.6). The
`build_and_solve(θ)` in `sensitivity_douglas_fir.jl` demonstrates the
state-dependent form for heat of sorption.

---

## 13.3 Local sensitivity by finite differences

A central difference of the whole solve approximates the Jacobian of an output
`y` with respect to `θ`:

```text
∂y/∂θ_k ≈ [ y(θ + h_k e_k) − y(θ − h_k e_k) ] / (2 h_k)
```

`FiniteDiff.finite_difference_jacobian` does the bookkeeping. With the
`p_inject` of §13.2:

```julia
using Pyrolysis, FiniteDiff

θ0 = [13.0, 1.88e5, 8.70e5]          # [log10(A / s⁻¹), E (J/mol), h (J/kg)]
base = p_inject(nothing, θ0)          # this p_inject ignores `base`

solve_opts = (radiation_model = SURFACE_ABSORPTION, abstol = 1e-9, reltol = 1e-7)
qoi(sol) = [sol.extras.total_mass[end],                     # mass left [kg]
            sum(sol.surface_T) / length(sol.surface_T)]     # mean surface T [K]
y(θ) = qoi(solve(p_inject(base, θ); progress = false, solve_opts...))

y0 = y(θ0)
J  = FiniteDiff.finite_difference_jacobian(y, θ0, Val(:central); relstep = 1e-4)
# J[i, k] = ∂y_i/∂θ_k; column 1 is per decade of A
```

- `y` must return a vector of the same length at every `θ`. Values on the fixed
  `output_times` grid are fine; anything whose length depends on the trajectory
  (the number of steps taken, a profile filtered of the `NaN` slots that
  depletion parking leaves) breaks the Jacobian assembly.
- For a scalar output use
  `FiniteDiff.finite_difference_gradient(f, θ, Val(:central); relstep)`.
- The step for entry `k` is `max(relstep * abs(θ_k), absstep)`, and `absstep`
  defaults to `relstep`. Pass `relstep` explicitly rather than relying on
  `FiniteDiff`'s default (about `6e-6` for central differences); §13.4 explains
  how to pick it.
- The cost is `2p` perturbed solves for `p = length(θ)`, plus the base solve.
  They are independent, so they parallelize (§13.5).

Pass every option that changes the solution — `radiation_model`, `use_ale`,
`jacobian`, tolerances — identically to every solve, here through
`solve_opts`.

---

## 13.4 Choosing the step

A finite difference balances two errors: truncation error, which grows as `h²`
for central differences, and the noise `ε` the time integration leaves in `y`,
which contributes about `ε/h`. On a well-posed problem the noise is small and
changes smoothly with `θ` — a small enough perturbation leaves the integrator's
step sequence unchanged — so the estimate settles over a wide range of steps.
On the §13.2 problem at `reltol = 1e-7`, the derivatives differ from those at
`relstep = 1e-4` by:

| Output | `relstep = 1e-2` | `1e-3` | `1e-5` to `1e-7` |
|---|---|---|---|
| peak surface temperature | 6 % | 0.06 % | ≤ 0.01 % |
| `T_s` at `t = 150 s` | 29 % | 0.4 % | ≤ 0.004 % |
| mean surface temperature | 0.2 % | 0.002 % | ≤ 0.0003 % |
| final mass, `total_mass[end]` | 0.06 % | 0.0006 % | ≤ 0.00002 % |

The large steps are truncation error, worst for values at one time, which move
with the reaction front. So:

1. **Start from `relstep = 1e-4`.**
2. **Check** by repeating with a step ten times larger and ten times smaller. If
   the estimates agree, keep them. If they scatter, the output is not smooth in
   `θ` at this point; look for one of the causes in §13.9.
3. **Do not ignore validation warnings.** A configuration `solve` warns about is
   often non-smooth in `θ` as well. The §13.2 material with its solids left at
   the default `λ = 0` (§13.9) needs 38 failed Newton iterations per solve, a
   `1e-6` relative change in `A` moves `T_s` by 0.35 K, and its derivatives
   change sign from one step to the next; with `λ = 1e-6` all of that is gone.

---

## 13.5 Running the perturbed solves in parallel

The `2p` perturbed solves are independent. `parameter_sweep` (Chapter 12) runs a
list of parameter vectors with `p_inject` as its `modify` and distributes them
over `Distributed` workers when any are present:

```julia
h   = max.(1e-4 .* abs.(θ0), 1e-4)                           # FiniteDiff's step, relstep 1e-4
pts = [θ0 .+ s .* h[k] .* (eachindex(θ0) .== k) for k in eachindex(θ0) for s in (1, -1)]
ys  = parameter_sweep(base, pts; modify = p_inject, extract = qoi, solve_opts...)
J   = reduce(hcat, [(ys[2k-1] .- ys[2k]) ./ (2h[k]) for k in eachindex(θ0)])
```

This builds the same central-difference Jacobian as §13.3. Define `p_inject`
and `qoi` with `@everywhere` when you run on workers (Chapter 12).

---

## 13.6 What the derivative means

A finite difference of `solve` is the derivative of the *discretized* model —
mesh, tolerances and all — not of the continuum equations. It changes with the
mesh like any other output, so take sensitivities on the mesh you calibrate on,
and recheck them if you refine it. Near a point where the output is not smooth
in `θ` (§13.9) the two one-sided differences disagree and the central difference
averages them; compare `y(θ + h e_k) − y(θ)` with `y(θ) − y(θ − h e_k)` when you
suspect one.

---

## 13.7 Inverse analysis: fitting parameters to data

Sensitivities support **inverse analysis** — recovering material parameters from
TGA, DSC, MCC, or cone-calorimeter measurements. The pattern is:

1. **Write `p_inject(base, θ)`** that rebuilds the problem with the candidate
   parameters at the right places (§13.2). Make sure `output_times` matches the
   times at which you have data so simulated and measured curves align.
2. **Write a scalar loss** comparing the relevant simulated field
   (`sol.surface_T`, `sol.mass_loss_rate`, `sol.extras.total_mass`, …) to data.
3. **Minimize it.** A derivative-free optimizer needs only the loss; a
   gradient-based one also needs a finite-difference gradient (`2p` solves per
   gradient).

The shipped fitting workflow, `examples/inverse_analysis/inverse_white_pine.jl`,
is derivative-free: `Optim.jl`'s `Fminbox(NelderMead())` on the loss, with a
`FiniteDiff` Hessian at the optimum for the standard errors (§14.7). It is robust
to the non-smooth responses of §13.9 and is the recommended starting point.

A gradient-based alternative hands `Optim` a central-difference gradient:

```julia
using Pyrolysis, FiniteDiff, Optim

T_data = load_my_surface_temperature()   # Vector{Float64} at output_times, your data
loss(θ) = sum(abs2, solve(p_inject(base, θ); progress = false,
                          solve_opts...).surface_T .- T_data)
g!(G, θ) = FiniteDiff.finite_difference_gradient!(G, loss, θ, Val(:central);
                                                  relstep = 1e-4)

θ0  = [13.0, 1.88e5, 8.70e5]             # [log10(A), E, h]
res = optimize(loss, g!, θ0, LBFGS())
θ_fit = Optim.minimizer(res)
```

Gradient methods stall when the loss is not smooth in `θ` (§13.9); switch to
`NelderMead()` if the line search keeps failing.

> **Scale your parameters.** Optimize `log10(A)` rather than `A`, and
> non-dimensionalize the rest if their magnitudes differ widely, so the optimizer
> sees comparable scales; otherwise the line search struggles and the gradient is
> lopsided.

The shipped `examples/df_sensitivity_comparison.jl` is **not** a programmatic
fitter — it parses a ThermaKin2Ds reference output
(`docs/DF_sensitivity_Cone_output.txt`) and an already-fitted Julia result CSV and
overlays them (MLR, surface/back temperature, mass, thickness, component MLR). Use
it as a template for *validating* a fit against a reference code, not for
generating one. The actual fitting workflows live under
`examples/inverse_analysis/`.

---

## 13.8 Global (Sobol) sensitivity via the sweep harness

Variance-based global sensitivity quantifies how much of an output's variance each
input (or group of inputs) is responsible for, over a *distribution* of inputs
rather than at a single point. Pyrolysis.jl has no dedicated Sobol function;
instead you build a sampling design, run it through `parameter_sweep`
(Chapter 12), and post-process. The shipped example
`examples/sensitivity_analysis/sensitivity_douglas_fir.jl` is a complete,
production-grade realization of this for a 16-parameter Douglas-fir cone model,
and `examples/sensitivity_analysis/sensitivity_sobol.jl` computes the variance
decomposition from its output.

The workflow has three stages:

**1. Sample the parameter space.** The example uses Owen-scrambled Sobol
quasi-random samples from `QuasiMonteCarlo.jl`, transformed to correlated
multivariate-Normal / Student-t groups whose covariances come from the
measurement campaign:

```julia
using QuasiMonteCarlo, Distributions, LinearAlgebra

# n_sobol unit-hypercube points, then map to each group's distribution.
sobol_unit = QuasiMonteCarlo.sample(
    n_samples, zeros(n_sobol), ones(n_sobol),
    SobolSample(R = OwenScramble(base = 2, pad = 32)),
)
# ... Φ⁻¹ to standard normal, Cholesky of the group covariance, t-scaling ...
# samples[:, j] is the j-th 16-parameter realization θ_j
```

**2. Solve every sample.** Each column of `samples` is a full parameter vector
`θ`; a `build_and_solve(θ)` function (the example's `p_inject`-equivalent)
rebuilds the Douglas-fir material and problem and solves it. Because every solve
is independent, this parallelizes cleanly with `Distributed.jl`:

```julia
using Distributed
@everywhere using Pyrolysis

# build_and_solve(θ) rebuilds the material/problem from θ and returns a solution.
# extract_qois(sol) reduces it to scalar QoIs + curves.
results = pmap(j -> extract_qois(build_and_solve(samples[:, j])), 1:n_samples)
```

This is the same separation as `parameter_sweep`: an immutable base, a `modify`/
`build_and_solve` that materializes a fresh problem per point, and an `extract`
that pulls out quantities of interest (peak MLR, time to ignition, total mass
loss, peak surface temperature, MLR/temperature curves, …). The example's solve
uses `use_ale = true`, `radiation_model = SURFACE_ABSORPTION`, and the default
Jacobian, which under ALE is `Structured(strategy = ApproxSolve())`.

**3. Decompose the variance.** `sensitivity_sobol.jl` fits a Polynomial Chaos
Expansion (PCE) surrogate to each scalar QoI and reads first-order and total
*group* Sobol indices analytically from the PCE coefficients, with bootstrap
confidence intervals. The Douglas-fir study groups the 16 parameters into 7
measurement sources (HFM/DSC calibration biases, three heat-capacity regressions,
heat-of-sorption parameters, and five heats of reaction), and reports how much of
each QoI's variance each *source* explains. A lighter-weight alternative used in
the same example is **Standardized Rank Regression Coefficients (SRRC)**, which
need only the sampled `(θ, QoI)` pairs and no surrogate.

> **You do not need PCE to get useful global sensitivity.** SRRC (rank regression
> of standardized inputs against the QoI) is a few lines on top of the sweep
> output and is robust for monotone responses. Use PCE when you want variance
> *decomposition* with interaction terms. Either way, the expensive part is the
> `N_samples` solves, which the harness parallelizes.

Global sampling spans parameter ranges far wider than a finite-difference step,
so the solver noise of §13.4 is small next to the response it measures.

For the theory of variance-based sensitivity and its relation to the sweep
harness, see Technical Reference §17.5.

---

## 13.9 Pitfalls

Sensitivities are only meaningful where the map `θ -> output` is well-behaved.
The package uses smooth constructions where it can (tanh temperature gates and
depletion limiter, smooth flux limiters — Technical Reference §17.6), but several
configurations still spoil a derivative:

- **Degenerate gas transport.** Condensed components default to `λ = 0`. Under
  Fickian gas transport, gas then cannot enter a cell that holds none: the gas
  field ends in a sharp front that advances only when integration noise seeds
  the next cell. The implicit Newton fails repeatedly on that front, the
  trajectory carries errors far above the tolerance whose timing changes from
  one `θ` to the next, and finite differences of it are meaningless (§13.4).
  `solve` warns in this configuration; give the condensed components a small
  `λ` (e.g. `1e-6`), as §13.2 does.

- **Warm start.** Reusing a prior solution as the initial state
  (`warm_start = previous_solution`) makes the result depend on history that
  `p_inject` does not rebuild. For sensitivity studies, **start each solve from
  the problem's own initial conditions** (do not pass `warm_start`), so every
  perturbed solve integrates the same map. This is why the sweep/sensitivity
  examples build a fresh problem per point and do not warm-start.

- **Non-smooth lateral-shrinkage laws (WithChi mode).** The cross-section law
  `A(t) = A_0 · law(χ̄)` is differentiated with ForwardDiff inside the residual
  (Technical Reference §10.6). A `law` with a kink (e.g. a hard `min`/`max` or a
  piecewise definition with a slope discontinuity) makes the output non-smooth in
  any parameter that moves `χ̄` across the kink. Use smooth shrink/swell laws.

- **Park events are discrete.** With `handle_depletion = true` a park rewrites
  the state in place; the trajectory is smooth in the parameters between events
  and jumps across one. A finite difference whose perturbed solves park at
  different times measures the jump. Keep depletion off in sensitivity work, or
  restrict the horizon to one without park events.

- **Tanh-gate leakage.** The smooth temperature gates leak a few percent of the
  reaction rate within roughly ±10 K of `T_min`/`T_max` (Technical Reference
  §6.2). Near those edges, `∂r/∂T_min` and `∂r/∂T_max` are smooth but small and
  not physically meaningful at the very edge of the window; interpret
  gate-boundary sensitivities with care.

---

## 13.10 Quick reference

| Task | Call |
|---|---|
| Jacobian of an output | `FiniteDiff.finite_difference_jacobian(θ -> y(solve(p_inject(base, θ))), θ0, Val(:central); relstep)` |
| Gradient of a scalar loss | `FiniteDiff.finite_difference_gradient(loss, θ0, Val(:central); relstep)` |
| Perturbed solves in parallel | `parameter_sweep(base, pts; modify = p_inject, extract)` (§13.5) |
| Global/Sobol | sample → `pmap`/`parameter_sweep` → PCE or SRRC (see examples) |

**Defaults to remember.** Start from `relstep = 1e-4` and check it with steps
ten times larger and smaller (§13.4); pass `log10(A)`, not `A`; heed `solve`'s
validation warnings (§13.9). For sensitivity-safe
solves: keep `warm_start = nothing`, `handle_depletion = false`, and the same
`solve` options for every perturbed solve.

For parameter sweeps (discrete response maps) see Chapter 12; for the solution
fields and accessors used in the outputs see Chapter 9; for the solver options
see Chapter 8; for the underlying theory see Technical Reference §17.
