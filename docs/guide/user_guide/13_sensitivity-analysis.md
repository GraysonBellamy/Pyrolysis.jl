# 13. Sensitivity Analysis

*Sensitivity analysis* asks how the outputs of a pyrolysis simulation change when
its inputs change. Pyrolysis.jl supports three complementary flavors: **local
forward** sensitivities (`∂output/∂θ` by automatic differentiation through the
whole solve), **local adjoint** sensitivities (`∂loss/∂θ` for a scalar loss, the
efficient mode for inverse problems with many parameters), and **global /
variance-based** (Sobol) sensitivities built on the parameter-sweep harness from
Chapter 12. This chapter shows how to use each, what can legitimately be a
parameter, how to wire your parameters into a fresh problem, and the numerical
pitfalls that make a sensitivity ill-defined. For the mathematics behind the
methods — Dual-number propagation, the augmented adjoint ODE, and the variance
decomposition — see Technical Reference §17.

> **Coordinate and symbol conventions.** `z = 0` is the bottom/substrate
> (boundary id 2, tag `:bottom`); `z = L` is the top/exposed surface (boundary id
> 1, tag `:top`); heat enters from the top. The Arrhenius parameters are the
> pre-exponential `A_i` [s⁻¹] (or [m³·kg⁻¹·s⁻¹] for higher-order reactions), the
> activation energy `E_i` [J·mol⁻¹], and the heat of reaction `h` [J·kg⁻¹] with
> the **`h > 0` endothermic** storage convention (Technical Reference §2 and §6).
> Mass concentration `ξ_j` [kg·m⁻³] is the primary species state.

---

## 13.1 Overview: when to use which method

| Method | Entry point | Returns | Cost scaling | Best when |
|---|---|---|---|---|
| **Local forward** | `forward_sensitivity` | output value + Jacobian `∂out/∂θ` | linear in `length(θ)` | few parameters (`≲ 20`), one or many outputs |
| **Local adjoint** | `adjoint_sensitivity` | scalar loss + gradient `∂loss/∂θ` | ~independent of `length(θ)` | many parameters, **one scalar** loss (calibration, optimization) |
| **Global (Sobol)** | `parameter_sweep` + a Sobol design | variance indices per parameter/group | `O(N_samples)` solves | uncertainty quantification, importance ranking over wide ranges |

Forward and adjoint are **local**: they linearize about a single parameter point
`θ`. Sobol is **global**: it explores a distribution of inputs and attributes
output variance to inputs, but it gives no derivative.

`forward_sensitivity` and `adjoint_sensitivity` are exported by `Pyrolysis` as
**stubs** — they error with a clear message until you load the
`PyrolysisSciMLSensitivityExt` extension (§13.6). The Sobol path needs no
extension; it is pure `parameter_sweep` plus your choice of sampling package.

---

## 13.2 The `p_inject` pattern: turning a vector into a problem

Both local methods share one idea. You hand them a **base problem**, a
**parameter vector `θ`**, and a **closure** that rebuilds a fresh
`PyrolysisProblem` from `(base_problem, θ)`:

```julia
p_inject(base_problem, θ) -> PyrolysisProblem
```

The extension differentiates the chain `θ -> output_fn(solve(p_inject(base, θ)))`
end-to-end. The whole solver (state cache, geometry, the boundary-temperature
Newton solve, the property derivatives) is type-generic in the scalar element
type, so when ForwardDiff passes a `Dual` vector the sensitivities propagate
through every stage. You never rewrite the residual.

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
differentiate:

```julia
using Pyrolysis

# θ = [A, E, h]
function p_inject(base, θ)
    A, E, h = θ[1], θ[2], θ[3]

    virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20)
    char   = SolidComponent(:char;   ρ = 200.0,  c = 1200.0, k = 0.10)
    gas    = GaseousComponent(:gas;   M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)

    rxn = Reaction(:decomp, 1 => (2, 3);
                   A = A, E = E, h = h,        # ← the differentiated parameters
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

A few rules that keep `p_inject` differentiable and consistent:

- **Place `θ` entries directly into the constructor calls.** The element type of
  `θ` flows into the `Reaction` and `SolidComponent` fields. Do not round, clamp,
  or `Float64(...)` a `θ` entry — that would break AD by stripping the `Dual`.
- **Everything not in `θ` is fixed.** Read fixed values from `base` if it is
  convenient, or just hard-code them as above.
- **Constructors that *must* receive `Float64`** (e.g. `tspan`, `output_times`,
  mesh thickness, `n_cells`) should be plain `Float64`/`Int` literals, never
  functions of `θ`. You differentiate physical parameters, not the time grid or
  the mesh count.
- **`yields` must sum to 1.0** for mass conservation (Technical Reference §6.5);
  `Reaction(...)` validates this at construction. If you make a yield a parameter
  you must rebuild the complementary yields so they still sum to 1.

### What can be a parameter `θ`

Anything that enters a constructor as a number can be a parameter. Common, useful
handles:

| Parameter | Where it lives | Notes |
|---|---|---|
| `A_i`, `E_i` | `Reaction(...; A=, E=)` | the dominant kinetic sensitivities |
| `h` (heat of reaction) | `Reaction(...; h=)` | scalar, tuple `(a,b)` → linear, or a closure |
| `ρ_j`, `c_p`, `k_j` | `SolidComponent(...; ρ=, c=, k=)` | pass a number, or wrap a `θ`-dependent `CallableProperty` |
| `λ_j`, `κ_j` | component `λ=`, `κ=` | gas-transfer / permeability for Darcy mode |
| `ε_j`, `α_j` | component `ε=`, `α=` | emissivity / mass-absorption |
| BC parameters | `HeatFluxBC`, `ConvectiveBC(h=,T_inf=)`, `RadiativeFluxBC(flux=)`, … | external flux, `h_conv`, `T_∞` |
| mixing `β` | `MixingSpec(WEIGHTED, β)` | conductivity / transport blend weight |

For a temperature-dependent property whose coefficients are parameters, build the
property as a closure over `θ` so the `Dual` propagates:

```julia
# θ = [a, b]; c_p(T) = a + b*T as a differentiable callable
cp = CallableProperty(T -> θ[1] + θ[2] * T)
virgin = SolidComponent(:virgin; ρ = 1190.0, c = cp, k = 0.20)
```

`Reaction` accepts a scalar `h` (a constant), a tuple `h = (a, b)` (linear in `T`,
`a + b*T`), or a property object such as `StateDependentProperty((T, ξ) -> ...)`
for moisture-dependent heats of sorption (Technical Reference §6.6). The
`build_and_solve(θ)` in `sensitivity_douglas_fir.jl` demonstrates the
state-dependent form for heat of sorption.

---

## 13.3 Forward-mode local sensitivity

`forward_sensitivity` differentiates an output (scalar or vector) with respect to
`θ` using `ForwardDiff.jacobian`. Its signature (attached by the extension to the
exported stub) is:

```julia
forward_sensitivity(
    base_problem,
    θ::AbstractVector,
    p_inject::Function;
    output_fn::Function = identity,   # solution -> output (scalar or vector)
    solve_kwargs...,                  # forwarded to every internal solve()
) -> (output, ∂output_∂θ)
```

- `output` is `output_fn(solve(p_inject(base, θ)))` evaluated at the base `θ`.
- `∂output_∂θ` is the Jacobian: if `output` is a length-`m` vector and `θ` has
  length `n`, the result is `m × n`. If `output_fn` returns a scalar, wrap it in a
  1-vector (ForwardDiff differentiates vector-valued functions).

`output_fn` reduces a `PyrolysisSolution` to the quantity you care about. Use the
solution's documented fields and accessors (Chapter 9): `sol.surface_T`,
`sol.mass_loss_rate`, `sol.T`, `sol.ξ`, `sol.extras.total_mass`,
`sol.extras.surface_mass_flux`, etc.

> **Cost.** ForwardDiff carries one extra "tangent" per parameter, so the solve
> cost scales roughly linearly in `length(θ)`. Forward mode is the right choice
> for a handful of parameters and is robust because it differentiates the exact
> code path the solver runs. For dozens of parameters and a single scalar loss,
> prefer adjoint mode (§13.5).

### Worked example: sensitivity of peak surface temperature to `(A, E, h)`

```julia
using Pyrolysis
using SciMLSensitivity, Zygote   # loads PyrolysisSciMLSensitivityExt

θ0 = [1.0e13, 1.88e5, 8.70e5]     # [A (1/s), E (J/mol), h (J/kg)]

# Scalar QoI wrapped in a 1-vector so ForwardDiff differentiates it.
peak_Ts(sol) = [maximum(sol.surface_T)]

# `base` is unused inside p_inject here, but it must still be passed
# (and could supply fixed pieces if you wire it in).
base = p_inject(nothing, θ0)

output, J = forward_sensitivity(base, θ0, p_inject;
                                output_fn = peak_Ts,
                                radiation_model = SURFACE_ABSORPTION,
                                abstol = 1e-9, reltol = 1e-7)

peak     = output[1]              # peak surface temperature [K] at θ0
dpeak_dA = J[1, 1]               # ∂(peak T_s)/∂A   [K·s]
dpeak_dE = J[1, 2]               # ∂(peak T_s)/∂E   [K·mol/J]
dpeak_dh = J[1, 3]              # ∂(peak T_s)/∂h   [K·kg/J]
```

`solve_kwargs` (here `radiation_model`, `abstol`, `reltol`) are forwarded
verbatim to each internal `solve` call; `progress = false` is forced internally.

### Differentiating a whole curve

`output_fn` may return a vector — e.g. the surface-temperature history sampled at
the problem's `output_times`. The Jacobian then has one row per saved time:

```julia
Ts_history(sol) = copy(sol.surface_T)             # length = number of saved times

output, J = forward_sensitivity(base, θ0, p_inject;
                                output_fn = Ts_history,
                                radiation_model = SURFACE_ABSORPTION)
# J[i, k] = ∂T_s(t_i)/∂θ_k
```

> **Keep the output dimension fixed.** `output_fn` must return a vector of the
> same length for every `θ` ForwardDiff probes. Returning `sol.surface_T`
> directly is fine because `saveat = output_times` is fixed and independent of
> `θ`. Returning something whose length depends on the trajectory (e.g. "number
> of time steps the integrator took", or — in ALE mode — a node array whose
> length changes after a depletion merge) will break the Jacobian assembly.

---

## 13.4 Manual forward sensitivity (without the extension)

`forward_sensitivity` is a thin convenience wrapper. If you would rather not pull
in `SciMLSensitivity`/`Zygote` just for forward mode, you can call `ForwardDiff`
yourself on the same `θ -> output_fn(solve(p_inject(base, θ)))` chain:

```julia
using Pyrolysis
using ForwardDiff

g(θ) = [maximum(solve(p_inject(base, θ);
                      progress = false,
                      radiation_model = SURFACE_ABSORPTION).surface_T)]

output = g(θ0)
J = ForwardDiff.jacobian(g, θ0)
```

This is exactly what the extension does internally; the extension just packages it
and provides the matching adjoint method. Forward mode needs only `ForwardDiff`
(a hard dependency of Pyrolysis), so this manual route works in any environment.

---

## 13.5 Adjoint-mode local sensitivity

For a **scalar loss** and **many parameters**, reverse mode (adjoint) computes the
full gradient at a cost roughly independent of `length(θ)`. The entry point is:

```julia
adjoint_sensitivity(
    base_problem,
    θ::AbstractVector,
    p_inject::Function,
    loss_fn::Function;                # solution -> scalar ℝ
    sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()),
    solve_kwargs...,
) -> (loss, ∂loss_∂θ)
```

- `loss_fn` reduces a solution to a single real number (a misfit, an objective).
- `sensealg` is a SciMLSensitivity reverse-mode algorithm. The default is
  `InterpolatingAdjoint(autojacvec = ReverseDiffVJP())`. It is forwarded to
  `solve` as the first-class `sensealg` keyword (the only reason `solve` accepts a
  `sensealg` kwarg at all), so the ODE integrator augments the solve with adjoint
  dynamics. The gradient is taken with `Zygote.gradient`.

### Worked example: gradient of a TGA/cone misfit

A typical calibration loss is the sum of squared errors between the simulated and
measured surface-temperature (or MLR) histories on the saved time grid:

```julia
using Pyrolysis
using SciMLSensitivity, Zygote

# Measured surface temperature at the same output_times as the problem.
T_data = load_my_surface_temperature()   # Vector{Float64}, your data

loss_fn(sol) = sum(abs2, sol.surface_T .- T_data)

θ0 = [1.0e13, 1.88e5, 8.70e5]
base = p_inject(nothing, θ0)

loss, grad = adjoint_sensitivity(base, θ0, p_inject, loss_fn;
                                 radiation_model = SURFACE_ABSORPTION,
                                 abstol = 1e-9, reltol = 1e-7)
# grad[k] = ∂loss/∂θ_k  — feed straight into a gradient-based optimizer
```

`grad` is what an optimizer (e.g. `Optim.LBFGS`, `Optimisers.Adam`) consumes to
fit `(A, E, h, …)` to data. This is the inverse-analysis recipe in §13.7.

> **Try forward mode first for small problems.** Adjoint mode carries
> significant machinery (continuous adjoint ODE, VJP construction) and can be
> sensitive to `sensealg` choice on stiff Arrhenius transients. For `≲ 20`
> parameters, forward mode is usually simpler and just as fast in wall-clock
> terms. Reserve adjoint mode for genuinely high-dimensional fits.

### Choosing a `sensealg`

`sensealg` is forwarded straight to the underlying ODE solve. SciMLSensitivity
offers several reverse-mode algorithms; the default `InterpolatingAdjoint` stores
a continuous interpolation of the forward solution and is a good general choice.
If you hit instability or excessive memory use, consult the SciMLSensitivity
documentation for alternatives (e.g. `QuadratureAdjoint`, `GaussAdjoint`, or a
different `autojacvec`); pass them through the `sensealg` keyword unchanged.

---

## 13.6 Loading the extension

`forward_sensitivity` and `adjoint_sensitivity` are exported stubs. Without the
extension they raise a clear error telling you which package to load. The
extension `PyrolysisSciMLSensitivityExt` is triggered automatically when **both**
`SciMLSensitivity` and `Zygote` are present in your environment:

```julia
using Pyrolysis
using SciMLSensitivity     # both are required to trigger the extension
using Zygote
```

The extension is declared in `Project.toml` as
`PyrolysisSciMLSensitivityExt = ["SciMLSensitivity", "Zygote"]`, so both packages
must be installed and loaded. (Forward mode via `ForwardDiff`, §13.4, does not
need the extension.) Once loaded, the concrete methods attach to the exported
stubs and the calls above work.

---

## 13.7 Inverse analysis: fitting parameters to data

The forward/adjoint machinery exists to support **inverse analysis** — recovering
material parameters from TGA, DSC, MCC, or cone-calorimeter measurements. The
pattern is:

1. **Write `p_inject(base, θ)`** that rebuilds the problem with the candidate
   parameters at the right places (§13.2). Make sure `output_times` matches the
   times at which you have data so simulated and measured curves align.
2. **Write a scalar `loss_fn(sol)`** comparing the relevant simulated field
   (`sol.surface_T`, `sol.mass_loss_rate`, `sol.extras.total_mass`, …) to data.
3. **Get the gradient** with `adjoint_sensitivity` (many parameters) or assemble
   it from `forward_sensitivity` on a 1-vector output (few parameters).
4. **Drive an optimizer** with `(loss, grad)`.

A compact LBFGS loop using forward-mode gradients (no extension required):

```julia
using Pyrolysis, ForwardDiff, Optim

T_data = load_my_surface_temperature()
loss(θ) = sum(abs2, solve(p_inject(base, θ);
                          progress = false,
                          radiation_model = SURFACE_ABSORPTION).surface_T .- T_data)

θ0 = [1.0e13, 1.88e5, 8.70e5]
res = optimize(loss, θ0, LBFGS(); autodiff = :forward)   # ForwardDiff under the hood
θ_fit = Optim.minimizer(res)
```

For larger parameter vectors, supply an explicit adjoint gradient instead:

```julia
function g!(G, θ)
    _, grad = adjoint_sensitivity(base, θ, p_inject,
                                  sol -> sum(abs2, sol.surface_T .- T_data);
                                  radiation_model = SURFACE_ABSORPTION)
    G .= grad
end
res = optimize(loss, g!, θ0, LBFGS())
```

> **Scale your parameters.** `A ~ 10¹³`, `E ~ 10⁵`, and `h ~ 10⁶` differ by many
> orders of magnitude. Optimize in log-space for `A` (`θ = log10(A)`) or
> non-dimensionalize so the optimizer sees comparable scales; otherwise the line
> search struggles and gradients are numerically lopsided. Apply the inverse
> transform inside `p_inject` (e.g. `A = 10^θ[1]`).

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
`Structured(strategy = DirectSolve())` Jacobian.

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

For the theory of variance-based sensitivity and its relation to the sweep
harness, see Technical Reference §17.5.

---

## 13.9 Pitfalls

Sensitivities are only meaningful where the map `θ -> output` is well-behaved.
The package is deliberately built to be AD-clean (smooth tanh temperature gates
and depletion limiter, midpoint-rule enthalpy, smooth flux limiters — Technical
Reference §17.6), but several configurations can still spoil a gradient:

- **Warm start.** Reusing a prior solution as the initial state
  (`warm_start = previous_solution`) makes the result depend on history, which is
  *not* differentiated by `p_inject`. For sensitivity studies, **start each solve
  from the problem's own initial conditions** (do not pass `warm_start`), so the
  differentiated chain is self-contained. This is why the sweep/sensitivity
  examples build a fresh problem per point and do not warm-start.

- **Non-smooth lateral-shrinkage laws (WithChi mode).** The cross-section law
  `A(t) = A_0 · law(χ̄)` is differentiated with ForwardDiff inside the residual
  (Technical Reference §10.6). A `law` with a kink (e.g. a hard `min`/`max` or a
  piecewise definition with a slope discontinuity) injects a non-smooth point into
  the sensitivity. Use smooth shrink/swell laws when you intend to differentiate
  through them.

- **Park events are discrete.** With `handle_depletion = true` a park rewrites
  the state in place; the trajectory is smooth in the parameters between events
  and jumps across one. Keep depletion off in sensitivity work, or restrict the
  horizon to one without park events.

- **Transpiration BC is forbidden with the default convention.** Leave
  `use_transpiration_bc = false` (the default). The volumetric advection source
  `S_conv` already carries gas enthalpy; enabling the surface transpiration term
  double-counts it, and problem construction / `solve` rejects it outright. There
  is no sensitivity question that requires turning it on.

- **Tanh-gate leakage.** The smooth temperature gates leak a few percent of the
  reaction rate within roughly ±10 K of `T_min`/`T_max` (Technical Reference
  §6.2). Near those edges, `∂r/∂T_min` and `∂r/∂T_max` are smooth but small and
  not physically meaningful at the very edge of the window; interpret
  gate-boundary sensitivities with care.

- **Integrator tolerances set a noise floor.** Finite `abstol`/`reltol` introduce
  small, parameter-dependent integration error. For tight sensitivities (and for
  reproducible adjoints), tighten tolerances (e.g. `abstol = 1e-9`,
  `reltol = 1e-7`, as the Douglas-fir example does) so the sensitivity signal sits
  well above the tolerance noise.

---

## 13.10 Quick reference

| Task | Call |
|---|---|
| Forward Jacobian of an output | `forward_sensitivity(base, θ, p_inject; output_fn)` |
| Forward Jacobian, no extension | `ForwardDiff.jacobian(θ -> out(solve(p_inject(base, θ))), θ)` |
| Adjoint gradient of a scalar loss | `adjoint_sensitivity(base, θ, p_inject, loss_fn; sensealg)` |
| Load the extension | `using SciMLSensitivity, Zygote` |
| Global/Sobol | sample → `pmap`/`parameter_sweep` → PCE or SRRC (see examples) |

**Defaults to remember.** `forward_sensitivity`: `output_fn = identity`.
`adjoint_sensitivity`: `sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP())`.
For sensitivity-safe solves: keep `warm_start = nothing`,
`handle_depletion = false`, `use_transpiration_bc = false`, the default
`jacobian = nothing` (→ `Structured(strategy = DirectSolve())`), and tighten
`abstol`/`reltol`.

For parameter sweeps (discrete response maps) see Chapter 12; for the solution
fields and accessors used in `output_fn`/`loss_fn` see Chapter 9; for the solver
options forwarded through `solve_kwargs` see Chapter 8; for the underlying theory
see Technical Reference §17.
