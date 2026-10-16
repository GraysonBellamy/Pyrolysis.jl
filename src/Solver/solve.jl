"""
Main solve interface for PyrolysisProblem.
"""

using ..Problem:
    PyrolysisProblem,
    PyrolysisSolution,
    validate_problem,
    initialize_state,
    StateLayout,
    get_temperature,
    get_concentration,
    get_node_positions,
    compute_volumes_from_state,
    _allocate_extras!
using ..Residual: residual!, cross_section_area_from_state
using ..Properties:
    effective_emissivity,
    effective_absorptivity,
    effective_absorption,
    effective_conductivity
using ..Discretization:
    compute_boundary_heat_fluxes, solve_boundary_temperature, _bc_boundary_face_fluxes
using ..Materials: has_darcy_flow
using ..BoundaryConditions:
    get_radiative_flux_bc,
    extract_radiative_flux_from_bcs,
    has_radiative_flux_bc,
    DirichletThermalBC,
    get_thermal_bc,
    get_mass_bc
using LinearSolve: KLUFactorization, UMFPACKFactorization  # the sparse defaults
using SparseArrays: nnz  # For logging sparse Jacobian info

import CommonSolve: solve  # Extend CommonSolve.solve for our problem type

# ============================================================================
# Post-processing helpers (private)
# ============================================================================
#
# The saved-output extraction needs the same per-species gas mass flux at the
# boundary that the live RHS computes: the blowing models of the convective
# BCs read it, so the boundary-temperature solve here reproduces what the
# integrator saw during the run.

# Per-species outward gas mass flux at a boundary face [kg/(m²·s)], positive
# leaving the slab at either boundary (`tag` is the boundary orientation): the
# Jacobian's boundary-face kernel (film law plus, in Darcy-Fick mode, the
# pressure-BC venting, mirroring the residual's mass BC), with its
# into-the-domain sign flipped. `0.0 - F` keeps a zero flux at +0.0.
@inline function _surface_gas_flux(
    material::Material{NC},
    ξ_face::NTuple{NC,Float64},
    T_face::Float64,
    mass_bc,
    d_half::Float64,
    t::Float64,
    darcy_active::Bool,
    tag::Val = Val(:top),
) where {NC}
    F = _bc_boundary_face_fluxes(
        T_face,
        ξ_face,
        mass_bc,
        material,
        d_half,
        t,
        tag,
        darcy_active,
    )
    return map(F_k -> 0.0 - F_k, F)
end

# ============================================================================
# Radiation Model Validation
# ============================================================================

"""
	_validate_radiation_model(model::RadiationModel, material)

Validate radiation model selection and warn about ignored parameters.

Issues warnings when:
- SURFACE_ABSORPTION is used with components that have α values (will be ignored)
- BEER_LAMBERT is used without any α values (no absorption will occur)
"""
function _validate_radiation_model(model::RadiationModel, material)
    # Warn if SURFACE_ABSORPTION is used with non-zero absorption coefficients
    if model == SURFACE_ABSORPTION
        if has_absorption_coefficients(material.components)
            @warn """
               SURFACE_ABSORPTION radiation model selected, but material components have
               non-zero absorption coefficients (α). These will be IGNORED.

               With SURFACE_ABSORPTION, all incident radiation is deposited in the surface
               cell regardless of material optical properties.

               If you want in-depth absorption using the specified α values, use:
               	radiation_model = BEER_LAMBERT
               """ maxlog = 1
        end
    end

    # Warn if BEER_LAMBERT is used without absorption coefficients
    if model == BEER_LAMBERT
        if !has_absorption_coefficients(material.components)
            @warn """
               BEER_LAMBERT radiation model selected, but no components have absorption
               coefficients (α) set. All radiation will pass through without absorption.

               Set absorption coefficients in your material definition:
               	SolidComponent(:name, ρ=..., c=..., k=..., α=1.47)  # m²/kg

               Or use SURFACE_ABSORPTION for opaque materials:
               	radiation_model = SURFACE_ABSORPTION
               """ maxlog = 1
        end
    end

    return nothing
end

# ============================================================================
# Reject any unknown kwargs so typos don't silently pass.
# ============================================================================

function _reject_unknown_kwargs(extra_kwargs)
    if !isempty(extra_kwargs)
        unknown = collect(keys(extra_kwargs))
        error(
            "solve(): unknown keyword arguments " *
            join(map(k -> "`$k`", unknown), ", ") *
            ". See `solve` docstring for the supported kwargs.",
        )
    end
    return nothing
end

"""
	_validate_solve_kwargs(; extra_kwargs, radiation_model,
	                          min_thickness, depletion_threshold, min_cells,
	                          diagnostics_stride, verify_jacobian)

Catch every cross-kwarg compatibility error at solve() entry, before any
allocation. Individual numeric or symbol checks belong here too whenever they
can be evaluated from kwargs alone (no context).
"""
function _validate_solve_kwargs(;
    extra_kwargs,
    radiation_model,
    min_thickness,
    depletion_threshold,
    min_cells,
    diagnostics_stride,
    verify_jacobian,
)
    # Unknown kwarg rejection.
    _reject_unknown_kwargs(extra_kwargs)

    # Cross-kwarg combinations.
    if radiation_model == P1_QUASI_STEADY
        error(
            "solve(): P1_QUASI_STEADY radiation lives in " *
            "`Pyrolysis.Experimental.RadiationP1` and is not part of the " *
            "core unified residual.",
        )
    end

    # Single-kwarg shape checks that fail fast.
    min_thickness > 0 ||
        throw(ArgumentError("solve(): min_thickness must be > 0, got $(min_thickness)."))
    (0 < depletion_threshold < 1) || throw(
        ArgumentError(
            "solve(): depletion_threshold must lie in (0, 1), got " *
            "$(depletion_threshold).",
        ),
    )
    min_cells >= 1 ||
        throw(ArgumentError("solve(): min_cells must be ≥ 1, got $(min_cells)."))
    diagnostics_stride >= 1 || throw(
        ArgumentError(
            "solve(): diagnostics_stride must be ≥ 1, got $(diagnostics_stride).",
        ),
    )

    # verify_jacobian shape: NamedTuple with at least :sample_times.
    if verify_jacobian !== nothing
        if !(verify_jacobian isa NamedTuple)
            throw(
                ArgumentError(
                    "solve(): verify_jacobian must be a " *
                    "NamedTuple `(rtol=, atol=, sample_times=)`, got " *
                    "$(typeof(verify_jacobian)).",
                ),
            )
        end
        allowed = (:rtol, :atol, :sample_times)
        bad = filter(k -> !(k in allowed), keys(verify_jacobian))
        isempty(bad) || throw(
            ArgumentError(
                "solve(): verify_jacobian unknown keys $(collect(bad)); " *
                "supported keys are $(allowed).",
            ),
        )
    end

    return nothing
end

# ============================================================================
# Auto-default Jacobian selection.
# ============================================================================
#
# On a fixed mesh the default is `Structured(strategy = DirectSolve())`: exact
# assembly of every coupling and a sparse direct solve.
#
# Under ALE the default is `Structured(strategy = ApproxSolve())`, with or
# without depletion handling. The `:full` ALE pattern carries the prefix-sum
# couplings of every node to every cell below it, which makes the factorized W
# effectively dense. The reduced pattern drops those couplings and keeps W
# local; the modified Newton iteration it drives still contracts, and the two
# strategies follow the same trajectory, through surface cells that collapse
# without parking too (Technical Reference §15.10.2).

_default_jacobian_backend(use_ale::Bool) =
    Structured(strategy = use_ale ? ApproxSolve() : DirectSolve())

# ============================================================================
# verify_jacobian DiscreteCallback
# ============================================================================
#
# When `solve(...; verify_jacobian = (rtol=, atol=, sample_times=))` is
# given, install a callback that fires at the listed sample times and
# compares the active backend's Jacobian against a `DenseAD_ForwardDiff`
# reference. Discrepancies log via `@warn` (do not throw — a single
# transient outlier should not abort an otherwise-good production run).

function _build_verify_callback(
    spec_active::JacobianSpec,
    def,
    ws,
    verify_kwargs::NamedTuple,
)
    rtol = get(verify_kwargs, :rtol, 1e-6)
    atol = get(verify_kwargs, :atol, 1e-10)
    sample_times = Float64.(collect(get(verify_kwargs, :sample_times, [0.0])))
    # Build the dense reference once and reuse it across firings.
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)

    # Behind a wrapper, like every problem-typed function the integrator
    # reaches (see `ErasedRHS`).
    verify_at = FunctionWrapper{Nothing,Tuple{_StateVector,Float64}}(
        (u, t) -> begin
            report = verify_jacobian_at(
                spec_active,
                def,
                ws,
                u,
                t;
                reference_spec = ref_spec,
                rtol = rtol,
                atol = atol,
            )
            if !report.passed
                @warn "verify_jacobian: backend $(typeof(spec_active.backend)) " *
                      "deviates from DenseAD reference at t=$(t): " *
                      "max_abs_dev=$(report.max_abs_dev), max_rel_dev=$(report.max_rel_dev), " *
                      "worst_index=$(report.worst_index)"
            end
        end,
    )

    condition =
        (u, t, integrator) -> begin
            # Fire once we cross any sample time. Tolerance: half the
            # integrator's current dt.
            dt = max(integrator.dt, eps(Float64))
            any(abs(t - st) <= 0.5 * dt for st in sample_times)
        end

    affect! = (integrator) -> begin
        verify_at(integrator.u, integrator.t)
        derivative_discontinuity!(integrator, false)
        return nothing
    end

    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# ============================================================================
# Linear solver / integrator selection
# ============================================================================

"""
	_resolve_integrator(; integrator, linear_solver, geometry_block=false)

Resolve the user's integrator / linear-solver request against the sparse
Jacobian prototype every `solve` attaches.

The default linear solver follows the Jacobian pattern: `UMFPACKFactorization`
when the pattern holds the ALE geometry block (`geometry_block`, see
[`materializes_geometry_block`](@ref)), and `KLUFactorization` for every other
pattern. UMFPACK's supernodal factorization is the fastest on the dense
geometry block. KLU refactors a fixed pattern in place without allocating and
is the fastest everywhere else, including the Beer–Lambert suffix
`DirectSolve` materializes on an Eulerian mesh, where UMFPACK's threshold
pivoting can return an inaccurate solve. The rule depends on the mode and
profile, not on the mesh size, so a configuration uses the same factorization
at every resolution.
"""
function _resolve_integrator(; integrator, linear_solver, geometry_block::Bool = false)
    default_ls = geometry_block ? UMFPACKFactorization() : KLUFactorization()
    ls = linear_solver === nothing ? default_ls : linear_solver
    if integrator === nothing
        return KenCarp4(autodiff = AutoFiniteDiff(), linsolve = ls)
    elseif integrator isa Type || integrator isa Function
        return integrator(autodiff = AutoFiniteDiff(), linsolve = ls)
    elseif linear_solver === nothing
        # Fully constructed algorithm.
        return integrator
    else
        error(
            "solve(): `linear_solver` can only be paired with `integrator` " *
            "when `integrator` is a constructor such as `KenCarp4`. To use " *
            "a fully constructed algorithm, pass it as `integrator = " *
            "KenCarp4(linsolve = ...)` on its own.",
        )
    end
end

# ============================================================================
# Main Solve Function
# ============================================================================

"""
	solve(problem::PyrolysisProblem; kwargs...) -> PyrolysisSolution

Solve a pyrolysis problem using DifferentialEquations.jl.

# Arguments
- `problem`: The PyrolysisProblem to solve. `solve` does not modify it: the run
  integrates a copy with its own mesh, so every solve of a problem, and every
  parameter-sweep point that reuses its mesh, starts from the same initial
  geometry. `solution.problem` is that copy; its mesh holds the run's final
  geometry, live-cell count and boundary faces.

# Keyword Arguments
- `integrator=nothing`: ODE integrator constructor (e.g. `KenCarp4`) or fully
  constructed algorithm (e.g. `KenCarp4(linsolve=KLUFactorization())`).
  `nothing` chooses `KenCarp4` with the default `linear_solver`.
- `linear_solver=nothing`: LinearSolve backend. `nothing` chooses by the
  Jacobian pattern (see "Linear Solver"): `UMFPACKFactorization()` for the
  exact ALE Jacobian, `KLUFactorization()` for every other Jacobian. Pair
  with `integrator` only when `integrator` is a constructor.
- `jacobian=nothing`: Jacobian backend. `nothing` picks
  `Structured(strategy = DirectSolve())` on a fixed mesh and
  `Structured(strategy = ApproxSolve())` under ALE. A `Structured` backend
  such as `Structured(strategy = DirectSolve())` overrides the default; the
  spec is built for this solve's own mesh copy, mode and options.
- `verify_jacobian=nothing`: NamedTuple `(rtol=, atol=, sample_times=)`
  installs a `DiscreteCallback` that compares the active backend against
  a `DenseAD_ForwardDiff` reference at the listed sample times. Reports
  via `@warn`; does not throw.
- `abstol=1e-8`: Absolute tolerance
- `reltol=1e-6`: Relative tolerance
- `sensealg=nothing`: Optional SciMLSensitivity sense algorithm forwarded to
  `OrdinaryDiffEq.solve`. This is primarily used by
  `PyrolysisSciMLSensitivityExt`.
- `maxiters=1000000`: Maximum iterations
- `dtmin=1e-12`: Minimum timestep
- `dtmax=Inf`: Maximum timestep
- `saveat=problem.output_times`: Times to save output
- `callback=nothing`: Additional callbacks
- `progress=true`: Show progress bar during solve
- `validate=true`: Validate problem before solving
- `use_ale=false`: Enable ALE mesh motion for shrinking/swelling materials
- `min_thickness=1e-6`: Minimum live thickness before termination (ALE only) [m]
- `handle_depletion=false`: Park depleted cells (see "Cell depletion" below)
- `depletion_threshold=0.05`: Depletion trigger. ALE: a cell is parked when
  its thickness falls below this fraction of its own initial thickness.
  Fixed mesh: the top live cell is parked when its condensed concentration
  falls below this fraction of its initial value.
- `min_cells=2`: Minimum number of live cells; the solve terminates instead
  of parking below it
- `warm_start=nothing`: Reuse a prior `PyrolysisSolution`'s final state as the initial state
  (`u0`). Component count, mesh size, and ALE/χ modes must match, and the prior solve must
  not have parked cells. `tspan[1]` still defines the integrator's start time, so the
  warm-started solve treats the prior final state as fresh initial data.

# Jacobian backend

Pass `jacobian = backend` to override the auto-default; the spec is built for
this solve's mode. The backend is `Structured`; behavior is selected by its
strategy:

| Spec                                                   | When to use                              |
|--------------------------------------------------------|------------------------------------------|
| `Structured(strategy = DirectSolve())`                 | Fixed-mesh default; exact assembly + sparse direct |
| `Structured(strategy = ApproxSolve())`                 | ALE default; drops the structured couplings (reduced pattern) |

# Returns
A `PyrolysisSolution` containing:
- Temperature and concentration profiles at saved times
- Mass loss rate history
- Surface temperature history
- Thickness history (when use_ale=true)
- Solver return code and timing

# ALE Mode
When `use_ale=true`, the solver uses the unified `ProblemDef`/`Workspace`
residual. Node positions evolve simultaneously with temperature and
concentration fields, maintaining proper conservation (GCL compliance). The
thickness history is extracted from the state vector and stored in
`solution.extras.thickness_history`.

# Cell Depletion Handling
When `handle_depletion=true`, in either mesh mode, the solver monitors
every live cell. Under ALE a cell that has shrunk below `depletion_threshold`
of its **own** initial thickness (a per-cell reference, so stretched meshes are
judged cell-by-cell) is parked; on a fixed mesh the top live cell is parked
once its condensed concentration has fallen below the same fraction of its
initial value. A fixed-mesh solve stops with a warning (retcode `:Terminated`)
when a cell below the top live cell becomes a void: the material has gasified
from within, the model carries no heat across the void, and a regressing
material belongs on the ALE mesh. A park merges the cell's mass and sensible enthalpy into its
neighbour and moves the freed state slot above the live surface as an inert
ghost: the state vector, mesh arrays and Jacobian keep their size, every
kernel runs over the live range `1:mesh.n_active_cells`, and the boundary
conditions act at the live surface. Vanishing cells therefore never enter the
nonlinear system (their `1/Δz²` conditioning is what stalls an implicit
solve), and the structured Jacobian stays valid through burnout. Ghost slots
read as `NaN` in `solution.T`, `solution.ξ`
and `solution.z`; `solution.extras.live_cells` and
`solution.extras.park_times` record the events. Parameter sensitivities are
valid only on horizons without park events.

# Recommended Solvers
The default is `KenCarp4(autodiff = AutoFiniteDiff())` with the default
Jacobian (see "Default Jacobian") and the default linear solver (see "Linear
Solver"). This is the production baseline for stiff pyrolysis runs.

Other stiff ODE algorithms such as `KenCarp3`, `TRBDF2`, `Rodas5P`, and
`QNDF` can be passed through `integrator`. Rosenbrock and BDF methods may
reevaluate or refactor the Jacobian more often on these problems, so benchmark
against the default before changing the solver family. Note: under
OrdinaryDiffEq v7 only the default solver set is re-exported by the umbrella
package, so passing e.g. `QNDF`/`KenCarp3`/`TRBDF2` requires loading the
relevant sublibrary (`OrdinaryDiffEqBDF`, `OrdinaryDiffEqSDIRK`, …) yourself.

`ILUGMRESFactorization` is available for large-ALE experiments where direct
factorization cost dominates, but it is not the default because inexact linear
solves can increase Newton rejections on strongly Arrhenius transients.

# Default Jacobian
With `jacobian = nothing` (the default), a fixed mesh uses
`Structured(strategy = DirectSolve())` and ALE uses
`Structured(strategy = ApproxSolve())`. The `Structured` backend assembles
per-operator contributions into a `StructuredJacobian`. `DirectSolve`
densifies the structured couplings (ALE prefix-sum, χ rank-1, Beer–Lambert
suffix-sum) at solve time, so the linear solve sees the exact Jacobian.
`ApproxSolve` drops them: under ALE the exact pattern couples every node to
every cell below it, and the reduced one keeps the Newton matrix local while
the modified Newton iteration still converges. Pass
`jacobian = Structured(strategy = DirectSolve())` for the exact ALE Jacobian.

# Linear Solver
With `linear_solver = nothing` the sparse factorization follows the Jacobian
pattern. `Structured(strategy = DirectSolve())` under ALE materializes the
geometry block, in which every cell couples to every node below it; that
pattern uses `UMFPACKFactorization()`, whose supernodal factorization handles
the dense block. Every other sparse Jacobian — every `ApproxSolve` Jacobian,
and `DirectSolve` on an Eulerian mesh, Beer–Lambert suffix included — uses
`KLUFactorization()`, which refactors the fixed pattern in place without
allocating. Pass `linear_solver = SparspakFactorization()`, another
LinearSolve backend, or a fully constructed `integrator` to choose the
factorization.

# Example
```julia
problem = PyrolysisProblem(mesh, material, bcs, (0.0, 100.0))
solution = solve(problem)
using OrdinaryDiffEqBDF  # QNDF lives in a sublibrary under OrdinaryDiffEq v7
solution = solve(problem, integrator=QNDF, abstol=1e-10)
solution = solve(problem, progress=false)  # Disable progress bar

# With ALE mesh motion for shrinking materials:
solution = solve(problem, use_ale=true, min_thickness=0.1e-3)
thickness_history = solution.extras.thickness_history

# Choose the factorization instead of the pattern-based default:
solution = solve(problem, use_ale=true, linear_solver=SparspakFactorization())

# The exact ALE Jacobian instead of the reduced default (factored by UMFPACK):
solution = solve(problem, use_ale=true, jacobian=Structured(strategy=DirectSolve()))

# With depletion handling (cells park at 1% of their initial thickness):
solution = solve(problem, use_ale=true, handle_depletion=true, depletion_threshold=0.01)

# Differential-correctness verification at sample times:
solution = solve(problem,
    verify_jacobian=(rtol=1e-6, atol=1e-10, sample_times=[0.0, 1.0, 10.0]))
```
"""
function solve(
    problem::PyrolysisProblem;
    linear_solver = nothing,
    integrator = nothing,

    # Pluggable Jacobian + verification:
    # `jacobian`: `nothing` → the mode's default (`DirectSolve` on a fixed
    # mesh, `ApproxSolve` under ALE); a `Structured` backend overrides it.
    # Either way the spec is built for this solve's (def, ws).
    jacobian::Union{Structured,Nothing} = nothing,
    # `verify_jacobian` installs a `DiscreteCallback` comparing the active
    # backend against a `DenseAD_ForwardDiff` reference at the listed
    # sample times. NamedTuple `(rtol=, atol=, sample_times=)` or `nothing`.
    verify_jacobian::Union{Nothing,NamedTuple} = nothing,
    # ---------------------------------------------------------------------

    abstol::Real = 1e-8,
    reltol::Real = 1e-6,
    sensealg = nothing,
    maxiters::Int = 1000000,
    dtmin::Real = 1e-12,
    dtmax::Real = Inf,
    saveat = problem.output_times,
    callback = nothing,
    progress::Bool = true,
    validate::Bool = true,
    use_ale::Bool = false,
    min_thickness::Real = 1e-6,
    # Surface cell depletion handling
    handle_depletion::Bool = false,
    depletion_threshold::Real = 0.05,
    min_cells::Int = 2,
    # Radiation model selection
    radiation_model::RadiationModel = NO_RADIATION,
    # Diagnostics log: a sensible-energy snapshot every `diagnostics_stride`
    # accepted steps, in `extras.diagnostics_log`.
    diagnostics::Bool = false,
    diagnostics_stride::Int = 10,

    # Reuse a prior `PyrolysisSolution`'s final state as the initial state.
    # Length and NC must agree with the new problem; `tspan[1]` still defines
    # the integrator's start time. Useful for parameter sweeps and continuation
    # studies where the previous solve already arrived at a near-converged state.
    warm_start::Union{Nothing,PyrolysisSolution} = nothing,

    # Catch any unknown kwargs via splat — `_reject_unknown_kwargs` raises a
    # clear error listing them.
    extra_kwargs...,
)
    _validate_solve_kwargs(;
        extra_kwargs = extra_kwargs,
        radiation_model = radiation_model,
        min_thickness = min_thickness,
        depletion_threshold = depletion_threshold,
        min_cells = min_cells,
        diagnostics_stride = diagnostics_stride,
        verify_jacobian = verify_jacobian,
    )

    # The caller's problem is left untouched; the run integrates a copy (see
    # `_problem_for_solve`), which the solution keeps as `solution.problem`.
    problem = _problem_for_solve(problem)

    # Validate problem
    if validate
        validate_problem(problem; strict = false)
    end

    # Get problem parameters (needed for validation)
    mesh = problem.mesh
    material = problem.material
    tspan = problem.tspan
    NC = n_components(material)

    # Validate radiation model selection
    _validate_radiation_model(radiation_model, material)

    # Start timing
    start_time = time()

    # Get radiation flux function from boundary conditions
    radiation_flux = _get_radiation_flux_function(problem)

    # Decide how the external radiation flux feeds the ODE.
    # - NO_RADIATION: no flux at all.
    # - SURFACE_ABSORPTION / BEER_LAMBERT: the flux is re-sampled at every RHS
    #   evaluation, so any time dependence reaches the volumetric deposit. One
    #   user-function call per RHS is negligible against the residual.
    uses_radiation =
        radiation_model == SURFACE_ABSORPTION || radiation_model == BEER_LAMBERT
    I_radiation = uses_radiation ? radiation_flux(0.0) : 0.0

    # Surface-transmissivity convention for the Beer-Lambert injection: on the
    # Kirchhoff default path (no explicit absorptivity on the RadiativeFluxBC)
    # the flux function returns the bare incident flux and the live
    # surface-cell ε_eff is applied at injection (α = ε, matching surface
    # emission). With an explicit absorptivity, the flux function has already
    # folded it in and injection must not rescale.
    rad_bc_top = get_radiative_flux_bc(get_thermal_bc(problem.bc_set, Val(:top)))
    radiation_kirchhoff = rad_bc_top === nothing || rad_bc_top.absorptivity === nothing

    # Under BEER_LAMBERT the volumetric sweep models only
    # top-incident radiation and the per-BC face flux is zeroed at every
    # boundary — a bottom RadiativeFluxBC would be silently dropped. Fail
    # loudly instead of losing the user's back-side heating.
    if radiation_model == BEER_LAMBERT &&
       has_radiative_flux_bc(problem.bc_set, Val(:bottom))
        throw(
            ArgumentError(
                "RadiativeFluxBC on the bottom boundary is not supported under " *
                "BEER_LAMBERT: the in-depth absorption sweep models top-incident " *
                "radiation only, and the boundary flux dispatch zeroes the BC, so " *
                "the back-side heating would be silently dropped. Use " *
                "SURFACE_ABSORPTION (face deposit) for back-side radiant heating, " *
                "or model it as a HeatFluxBC.",
            ),
        )
    end

    if radiation_model == SURFACE_ABSORPTION
        @info "Using SURFACE_ABSORPTION radiation model"
    elseif radiation_model == BEER_LAMBERT
        @info "Using BEER_LAMBERT radiation model"
    end

    # Canonical setup: ICs → ProblemDef → Workspace → u0 (single sequence,
    # shared with the exported `initialize_state`).
    def, ws, u0 = initialize_state(
        problem;
        use_ale = use_ale,
        radiation_model = radiation_model,
        radiation_intensity = I_radiation,
        radiation_kirchhoff = radiation_kirchhoff,
    )

    # Create solution container (after initialize_state so any mesh state it
    # snapshots reflects the applied initial conditions).
    solution = PyrolysisSolution(problem)

    if warm_start !== nothing
        _apply_warm_start!(u0, problem, ws.layout, warm_start)
    end
    rhs! =
        let def = def,
            ws = ws,
            radiation_flux = radiation_flux,
            uses_radiation = uses_radiation

            (du, u, p, t) -> begin
                if uses_radiation
                    ws.I_radiation = radiation_flux(t)
                end
                residual!(du, u, def, t, ws)
            end
        end

    # Build the Jacobian spec for this solve's own (def, ws).
    backend = jacobian === nothing ? _default_jacobian_backend(use_ale) : jacobian
    jac_spec = JacobianSpec(backend, def, ws)
    # The Jacobian refreshes the incident flux at its own `t`, as the RHS does.
    jac! =
        let s = jac_spec,
            def = def,
            ws = ws,
            radiation_flux = radiation_flux,
            uses_radiation = uses_radiation

            (J, u, p, t) -> begin
                if uses_radiation
                    ws.I_radiation = radiation_flux(t)
                end
                compute_jacobian!(J, s, u, def, t, ws)
            end
        end

    # The union of physics + geometry sparsity is the allocation template.
    # Both are `SparseMatrixCSC{Bool,Int}`; the Float64 working matrix is
    # materialized by OrdinaryDiffEq.
    S_geom = jac_spec.sparsity_geom
    S_union =
        S_geom === nothing ? jac_spec.sparsity_phys : (jac_spec.sparsity_phys .| S_geom)
    jac_prototype = Float64.(S_union)

    # The pattern fixes the default factorization (UMFPACK when it holds the
    # ALE geometry block, KLU otherwise), so the integrator is resolved after
    # the spec is built.
    geometry_block =
        materializes_geometry_block(def.mode, sparsity_profile(jac_spec.backend.strategy))
    solver_alg = _resolve_integrator(; integrator, linear_solver, geometry_block)

    # The right-hand side and the Jacobian reach the integrator behind
    # problem-independent types (`ErasedRHS`, `ErasedJacobian`), as do the
    # callbacks below: the integrator then compiles once per solver
    # configuration rather than once per material or boundary-condition type.
    # `FullSpecialize`: under the default `AutoSpecialize` the integrator stores its
    # callbacks type-erased, and inference then compiles the stepper for a second,
    # callback-free integrator type as well, which doubles its compile time.
    ode_fn_with_jac = OrdinaryDiffEq.ODEFunction{true,FullSpecialize}(
        ErasedRHS(rhs!);
        jac = ErasedJacobian{typeof(jac_prototype)}(jac!),
        jac_prototype = jac_prototype,
    )
    ode_problem = ODEProblem(ode_fn_with_jac, u0, tspan)

    linsolve_note =
        hasproperty(solver_alg, :linsolve) ?
        " and `$(nameof(typeof(solver_alg.linsolve)))`" : ""
    @info "Using `$(typeof(jac_spec.backend))` Jacobian " *
          "with sparse storage ($(nnz(jac_prototype)) non-zeros)$(linsolve_note)"
    if hasproperty(solver_alg, :linsolve) && solver_alg.linsolve isa ILUGMRESFactorization
        @warn "Using ILU-GMRES iterative solver (τ=$(solver_alg.linsolve.τ)). " *
              "Its inexact linear solves can add Newton failures and rejected " *
              "steps on strongly Arrhenius transients; the default direct " *
              "factorization solves each Newton system exactly." maxlog = 1
    end

    # Build solver options. The tolerances, step bounds and output times enter
    # the integrator's type, so they are converted to one type each.
    solve_kwargs = Dict{Symbol,Any}(
        :abstol => Float64(abstol),
        :reltol => Float64(reltol),
        :maxiters => maxiters,
        :dtmin => Float64(dtmin),
        :dtmax => Float64(dtmax),
        :dt => problem.dt_initial,
        :progress => false,  # Disable OrdinaryDiffEq's progress, we use our own
    )
    if sensealg !== nothing
        solve_kwargs[:sensealg] = sensealg
    end

    solve_kwargs[:saveat] = saveat isa Number ? Float64(saveat) : collect(Float64, saveat)

    # Build callback list
    callbacks_list = []

    # StepUpdateCallback — mirror integrator state into the mutable mesh at the
    # end of each accepted step. Residual / Jacobian assembly
    # is mesh-read-only, so we need an explicit sync point for AMR, depletion
    # checks, and post-processing readers. Push first so it runs before any
    # downstream callback that reads mesh.
    push!(callbacks_list, create_step_update_callback(ws, material))

    # Depletion bookkeeping shared by the depletion callback, the thickness
    # termination callback and solution extraction.
    depletion_cb_state =
        handle_depletion ? DepletionCallbackState(ws.layout, mesh, material, u0) : nothing

    # Thickness termination (ALE): the live thickness z[n_live+1] − z[1].
    if use_ale
        sv_ref = depletion_cb_state !== nothing ? depletion_cb_state : ws.layout
        thickness_cb = create_thickness_termination_callback(sv_ref, Float64(min_thickness))
        push!(callbacks_list, thickness_cb)
    end

    if depletion_cb_state !== nothing
        depletion_cb = create_depletion_callback(
            ws,
            mesh,
            material,
            depletion_cb_state;
            threshold = Float64(depletion_threshold),
            min_cells = min_cells,
        )
        push!(callbacks_list, depletion_cb)
    end

    # Newton-failure step cap for integrators that solve their stages by Newton.
    if hasproperty(solver_alg, :nlsolve)
        push!(callbacks_list, create_newton_failure_cap_callback())
    end

    # Progress bar, advanced after every accepted step. A callback that fires at
    # fixed times would force the integrator to stop at each of them. The
    # callback is installed with the bar disabled too, so both settings give
    # the integrator the same type and share its compiled code.
    t_start, t_end = tspan
    t_span_length = t_end - t_start
    last_progress = Ref(0)
    prog = Progress(100; desc = "Solving: ", showspeed = true, enabled = progress)
    progress_cb = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> begin
            current_progress = floor(Int, 100 * (integrator.t - t_start) / t_span_length)
            if current_progress > last_progress[]
                update!(prog, current_progress)
                last_progress[] = current_progress
            end
            derivative_discontinuity!(integrator, false)
        end;
        save_positions = (false, false),
    )
    push!(callbacks_list, progress_cb)

    # Add diagnostics callback if requested. Snapshot mode — records the
    # sensible energy at stride boundaries; the cumulative ledger fields of the
    # row stay NaN.
    diagnostics_log = nothing
    if diagnostics
        diagnostics_log = DiagnosticsRow[]
        stride_count = Ref(0)
        mesh_for_diag = ws.mesh
        material_for_diag = problem.material
        diagnostics_record! =
            let log = diagnostics_log,
                stride_count = stride_count,
                stride = diagnostics_stride,
                mesh = mesh_for_diag,
                material = material_for_diag,
                layout = ws.layout

                (u, t) -> begin
                    stride_count[] += 1
                    if stride_count[] % stride == 0
                        # Sync mesh state from the integrator's u so the snapshot
                        # reads the current cell states.
                        update_mesh_from_state!(mesh, u, layout)
                        push!(log, snapshot_diagnostics_row(t, mesh, material))
                    end
                    return nothing
                end
            end
        diagnostics_cb = DiagnosticsCallback(diagnostics_record!)
        push!(callbacks_list, diagnostics_cb)
    end

    # verify_jacobian DiscreteCallback: fires at the requested sample times
    # and warns (does not throw) if the active backend deviates from the
    # DenseAD reference beyond tolerance.
    if verify_jacobian !== nothing
        push!(callbacks_list, _build_verify_callback(jac_spec, def, ws, verify_jacobian))
    end

    # Add user callback if provided
    if callback !== nothing
        push!(callbacks_list, callback)
    end

    # Combine callbacks
    if !isempty(callbacks_list)
        if length(callbacks_list) == 1
            solve_kwargs[:callback] = callbacks_list[1]
        else
            solve_kwargs[:callback] = CallbackSet(callbacks_list...)
        end
    end

    # Use sparse solver if configured, otherwise user-provided or default

    # Solve
    raw_solution = OrdinaryDiffEq.solve(ode_problem, solver_alg; solve_kwargs...)

    progress && finish!(prog)

    # Record solve time
    solution.solve_time = time() - start_time

    # Hand the diagnostics log over to the solution (empty when not enabled).
    if diagnostics_log !== nothing
        append!(solution.extras.diagnostics_log, diagnostics_log)
    end

    # Extract and store results based on ALE mode
    if use_ale
        extract_solution_ale!(
            solution,
            raw_solution,
            problem,
            ws.layout;
            depletion = depletion_cb_state,
            radiation_model = radiation_model,
        )
    else
        extract_solution!(
            solution,
            raw_solution,
            problem;
            depletion = depletion_cb_state,
            radiation_model = radiation_model,
        )
    end

    return solution
end

# ============================================================================
# Helper Functions
# ============================================================================

# The problem a solve integrates: the caller's, with its own copy of the mesh.
# The initial conditions, the dry-solid reference masses, the per-step mesh
# sync, parks and the boundary faces they move all write that copy, so every
# solve starts from the mesh the caller built, however often it is solved.
_problem_for_solve(p::PyrolysisProblem) = PyrolysisProblem(
    deepcopy(p.mesh),
    p.material,
    p.bc_set,
    p.T_initial,
    p.ξ_initial,
    p.tspan,
    p.dt_initial,
    p.output_times,
)

"""
	_apply_warm_start!(u0, problem, layout, prev::PyrolysisSolution)

Overwrite `u0` with the final raw state from a prior `PyrolysisSolution`.

Compatibility checks:
- `prev.raw_solution` must exist (the prior solve must have populated it).
- The previous solve must not have parked cells: its final state holds ghost
  slots, and a new solve starts with every cell live.
- The previous problem's component count must match `n_components(layout)`.
- The previous final state vector length must match `length(u0)`. This implicitly
  enforces matching mesh size, ALE/Eulerian mode, and χ mode, since those drive
  state length together with NC.
"""
function _apply_warm_start!(
    u0::AbstractVector,
    problem::PyrolysisProblem,
    layout,
    prev::PyrolysisSolution,
)
    if prev.raw_solution === nothing
        error(
            "warm_start: previous solution has no raw_solution. " *
            "The prior solve() may have failed before extraction.",
        )
    end
    if !isempty(prev.extras.park_times)
        error(
            "warm_start: the previous solve parked $(length(prev.extras.park_times)) " *
            "cell(s), so its final state holds ghost slots, while a new solve starts " *
            "with every cell live. Warm-start only from a solve without park events.",
        )
    end
    prev_NC = n_components(prev.problem.material)
    cur_NC = n_components(problem.material)
    if prev_NC != cur_NC
        error(
            "warm_start: previous solution has $(prev_NC) components, " *
            "current problem has $(cur_NC). Cannot reuse state.",
        )
    end
    prev_u = last(prev.raw_solution.u)
    if length(prev_u) != length(u0)
        error(
            "warm_start: state vector length mismatch (previous = " *
            "$(length(prev_u)), current = $(length(u0))). Mesh size, ALE " *
            "mode, or χ mode must agree.",
        )
    end
    @inbounds for i in eachindex(u0)
        u0[i] = prev_u[i]
    end
    return u0
end

"""
	_get_radiation_flux_function(problem::PyrolysisProblem) -> Function

Get the radiation flux function from `RadiativeFluxBC` on the top boundary.

Returns a function t -> I [W/m²]: the bare incident flux on the Kirchhoff
default path (the live surface-cell ε_eff is applied at injection by
`compute_beer_lambert_source!`), or the flux with the BC's explicit
absorptivity already folded in (see `extract_radiative_flux_from_bcs`).
When no `RadiativeFluxBC` is present, returns `t -> 0.0` (no radiative heating).

The radiation_model passed to solve() determines how the flux is used:
- SURFACE_ABSORPTION: flux enters via the surface energy balance (boundary Newton)
- BEER_LAMBERT: flux used for volumetric Beer-Lambert absorption (I_radiation)
"""
function _get_radiation_flux_function(problem::PyrolysisProblem)
    if has_radiative_flux_bc(problem.bc_set, Val(:top))
        return function (t)
            I, _ = extract_radiative_flux_from_bcs(problem.bc_set, Val(:top), t)
            return I
        end
    else
        return _zero_radiation_flux
    end
end

@inline _zero_radiation_flux(_t) = 0.0

"""
	extract_solution!(solution::PyrolysisSolution, raw_solution, problem; radiation_model)

Extract results from raw ODE solution into PyrolysisSolution structure.
"""
function extract_solution!(
    solution::PyrolysisSolution,
    raw_solution,
    problem;
    depletion::Union{Nothing,DepletionCallbackState} = nothing,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
)
    mesh = problem.mesh
    material = problem.material
    NC = n_components(material)
    NR = n_reactions(material)
    n_cells_val = length(mesh.cells)

    # Store raw solution for interpolation
    solution.raw_solution = raw_solution

    # Get return code
    solution.retcode = Symbol(raw_solution.retcode)

    # Extract times
    times = raw_solution.t
    solution.t = copy(times)

    # Pre-allocate output arrays. Parked ghost slots read as NaN.
    n_times = length(times)
    solution.T = fill(NaN, n_cells_val, n_times)
    solution.ξ = fill(NaN, NC, n_cells_val, n_times)
    solution.z = Matrix{Float64}(undef, 0, 0)             # Eulerian: no node positions
    solution.surface_T = zeros(n_times)
    solution.mass_loss_rate = zeros(n_times)
    _allocate_extras!(solution, n_times)
    depletion === nothing || append!(solution.extras.park_times, depletion.park_times)

    # Variable cross-section diagnostics — non-ALE has no χ evolution, so
    # A(t) ≡ cross_section_area_initial and MLR_gauge ≡ MLR_local.
    @assert mesh.cross_section_area == mesh.cross_section_area_initial
    fill!(solution.extras.cross_section_area_history, mesh.cross_section_area_initial)

    # Pre-allocate work arrays for reaction rate computation
    rates = zeros(NR)
    heats = zeros(NR)

    # Get boundary conditions for flux computation
    bc_set = problem.bc_set
    thermal_bc_top = get_thermal_bc(bc_set, Val(:top))
    thermal_bc_bottom = get_thermal_bc(bc_set, Val(:bottom))
    mass_bc_top = get_mass_bc(bc_set, Val(:top))
    mass_bc_bottom = get_mass_bc(bc_set, Val(:bottom))

    # Extract data at each saved time
    for (i, t) in enumerate(times)
        u = raw_solution.u[i]
        n_live = depletion === nothing ? n_cells_val : n_live_cells(depletion, t)
        solution.extras.live_cells[i] = n_live
        solution.extras.thickness_history[i] =
            mesh.face_positions[mesh.cells[n_live].face_ids[2]] -
            mesh.face_positions[mesh.cells[1].face_ids[1]]

        # Extract temperatures of the live cells (ghost slots stay NaN)
        T_profile = zeros(n_live)
        @inbounds for j = 1:n_live
            T_profile[j] = u[j]
            solution.T[j, i] = T_profile[j]
        end

        # Extract concentrations
        ξ_profile = zeros(NC, n_live)
        @inbounds for j = 1:n_live
            for k = 1:NC
                ξ_profile[k, j] = u[n_cells_val+(k-1)*n_cells_val+j]
                solution.ξ[k, j, i] = ξ_profile[k, j]
            end
        end

        # ================================================================
        # Solve for true boundary temperatures (not cell-center values).
        # The per-species gas mass flux at each face is reconstructed the
        # same way the live mass BC application does, so the blowing models
        # see the flux the integrator saw.
        # ================================================================

        # Surface (top, boundary_id=1) — the live top cell
        T_cell_top = T_profile[n_live]
        ξ_top = ntuple(k -> ξ_profile[k, n_live], NC)
        k_top = effective_conductivity(ξ_top, T_cell_top, material)
        ε_top = effective_emissivity(ξ_top, T_cell_top, material)
        a_top = effective_absorptivity(ξ_top, T_cell_top, material)
        Δz_half_top =
            mesh.cell_volumes[n_live] / mesh.face_areas[mesh.cells[n_live].face_ids[2]] / 2

        J_top = _surface_gas_flux(
            material,
            ξ_top,
            T_cell_top,
            mass_bc_top,
            Δz_half_top,
            t,
            has_darcy_flow(material),
        )

        state_top = solve_boundary_temperature(
            T_cell_top,
            k_top,
            Δz_half_top,
            thermal_bc_top,
            material,
            J_top,
            t,
            ε_top;
            a_eff = a_top,
            radiation_model = radiation_model,
        )
        solution.surface_T[i] = state_top.T_surface
        solution.extras.surface_blowing_ratio[i] = state_top.blowing_ratio
        solution.extras.surface_pooling_factor[i] = state_top.pooling_factor
        solution.extras.surface_face_fraction[i] = state_top.face_fraction
        solution.extras.surface_face_buoyancy[i] = state_top.face_buoyancy

        # Back (bottom, boundary_id=2) — vents through the user-supplied bottom
        # mass BC (zero for the default impermeable holder; the RHS applies
        # a permeable back-face BC, so the MLR must count it).
        T_cell_bottom = T_profile[1]
        ξ_bottom = ntuple(k -> ξ_profile[k, 1], NC)
        k_bottom = effective_conductivity(ξ_bottom, T_cell_bottom, material)
        ε_bottom = effective_emissivity(ξ_bottom, T_cell_bottom, material)
        a_bottom = effective_absorptivity(ξ_bottom, T_cell_bottom, material)
        Δz_half_bottom =
            mesh.cell_volumes[1] / mesh.face_areas[mesh.cells[1].face_ids[1]] / 2
        J_bottom = _surface_gas_flux(
            material,
            ξ_bottom,
            T_cell_bottom,
            mass_bc_bottom,
            Δz_half_bottom,
            t,
            has_darcy_flow(material),
            Val(:bottom),
        )
        state_bottom = solve_boundary_temperature(
            T_cell_bottom,
            k_bottom,
            Δz_half_bottom,
            thermal_bc_bottom,
            material,
            J_bottom,
            t,
            ε_bottom;
            a_eff = a_bottom,
            radiation_model = radiation_model,
        )
        solution.extras.back_temperature[i] = state_bottom.T_surface

        # Total mass (integrate ξ over the live cells)
        total_m = 0.0
        for j = 1:n_live
            cell_vol = mesh.cell_volumes[j]
            for k = 1:NC
                total_m += ξ_profile[k, j] * cell_vol
            end
        end
        # Normalize by original area for cone-gauge-style mass loading [kg/m²]
        solution.extras.total_mass[i] = total_m / mesh.cross_section_area_initial

        # ================================================================
        # Boundary flux computation (using solved surface temperature)
        # ================================================================

        # Surface cell properties (top boundary, z=L, boundary_id=1)
        # Note: ξ_top, k_top, ε_top already computed above for boundary temperature solve
        T_surface_solved = state_top.T_surface

        # Store effective properties
        solution.extras.surface_emissivity[i] = ε_top
        solution.extras.surface_conductivity[i] = k_top

        # Compute heat flux components using SOLVED surface temperature.
        # The face the Newton applied the blowing models to, so the convective
        # decomposition matches the corrected flux.
        heat_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_top,
            T_surface_solved,
            t,
            ε_top,
            radiation_model,
            BlowingFace(J_top, T_surface_solved, material);
            a_eff = a_top,
        )
        solution.extras.surface_heat_flux_absorbed[i] = heat_fluxes.absorbed
        solution.extras.surface_heat_flux_reradiation[i] = heat_fluxes.reradiation
        solution.extras.surface_heat_flux_convection[i] = heat_fluxes.convection
        # Dirichlet: the applied face flux is the conducted flux from the
        # boundary solve (the BC itself has no closed-form flux).
        solution.extras.surface_heat_flux_net[i] =
            thermal_bc_top isa DirichletThermalBC ? state_top.q_conducted : heat_fluxes.net

        # Conductive heat flux at surface: use solved boundary state
        # This is consistent with the boundary temperature solve and ensures energy balance
        solution.extras.surface_heat_flux_conduction[i] = state_top.q_conducted

        # Volumetric radiative absorption (BEER_LAMBERT only): the in-depth
        # deposit Σ q'''·Δz telescopes to α_surf·q_ext·(1 − e^{−τ_tot}) —
        # the energy-ledger companion to the face-applied
        # `surface_heat_flux_net`, which excludes it. α_surf is the BC's
        # explicit absorptivity (already folded into the extracted flux) or
        # the surface cell's live ε_eff on the Kirchhoff default path.
        if radiation_model == BEER_LAMBERT
            I_inc, has_explicit_α = extract_radiative_flux_from_bcs(bc_set, Val(:top), t)
            I_eff = has_explicit_α ? I_inc : ε_top * I_inc
            τ_tot = 0.0
            for j = 1:n_live
                ξ_j = ntuple(k -> ξ_profile[k, j], NC)
                τ_tot += effective_absorption(ξ_j, material) * Geometry.cell_size(mesh, j)
            end
            solution.extras.volumetric_absorption[i] = I_eff * (1.0 - exp(-τ_tot))
            solution.extras.radiation_transmitted[i] = I_eff * exp(-τ_tot)
        end

        # Surface MLR: per-species gas flux already used in the boundary
        # solve. `surface_mass_flux` is the top face only; the total MLR also
        # counts the bottom film/venting flux (zero for the default
        # impermeable holder).
        total_mass_flux = sum(J_top)
        bottom_mass_flux = sum(J_bottom)
        solution.extras.surface_mass_flux[i] = total_mass_flux

        # Variable cross-section MLR diagnostics — non-ALE keeps A(t) ≡ A_0
        # so the gauge and local flavors coincide.
        solution.mass_loss_rate[i] = total_mass_flux + bottom_mass_flux
        solution.extras.MLR_gauge[i] = total_mass_flux + bottom_mass_flux

        # Back face heat flux (boundary_id=2, z=0) - use solved boundary temperature
        T_back_solved = state_bottom.T_surface
        back_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_bottom,
            T_back_solved,
            t,
            ε_bottom,
            radiation_model,
            BlowingFace(J_bottom, T_back_solved, material);
            a_eff = a_bottom,
        )
        solution.extras.back_heat_flux[i] =
            thermal_bc_bottom isa DirichletThermalBC ? state_bottom.q_conducted :
            back_fluxes.net

        # Gas generation rate: compute from reaction rates (physics-based)
        # This is NOT the actual mass loss rate (surface flux), but rather the
        # volumetric production of gases throughout the sample.
        # Total gas generation = Σ_cells Σ_reactions (gas_produced * r * V)
        # This is much smoother than numerical differentiation of total mass
        gen_rate = 0.0
        for j = 1:n_live
            T = T_profile[j]
            ξ = ntuple(k -> ξ_profile[k, j], NC)
            cell_vol = mesh.cell_volumes[j]

            # Evaluate reaction rates at this cell
            Physics.evaluate_reactions!(rates, heats, ξ, T, material)

            # Sum gas production from each reaction
            for r_idx = 1:NR
                rxn = material.reactions[r_idx]
                r = rates[r_idx]

                # Gas production: sum mass yields for any gas-phase product
                for (p_idx, ν) in rxn.products
                    if material.components[p_idx] isa GaseousComponent
                        gen_rate += ν * r * cell_vol
                    end
                end
            end
        end
        # Normalize by original area for cone-gauge-style generation rate [kg/(m²·s)]
        solution.extras.gas_generation_rate[i] = gen_rate / mesh.cross_section_area_initial
    end

    _warn_transmitted_radiation(solution, radiation_model)
    return solution
end

# Beer–Lambert radiation that reaches the back face leaves the domain: the
# sweep deposits nothing there and the back-face boundary condition sees only
# the surface it describes. That is the right answer for a slab whose back
# face really is open, and a modelling error for one that sits on an opaque
# holder with a transparent absorption coefficient, so the extracted
# `radiation_transmitted` is reported and a run that loses more than 1 % of
# the entering flux at any saved time says so once.
function _warn_transmitted_radiation(solution::PyrolysisSolution, radiation_model)
    radiation_model == BEER_LAMBERT || return nothing
    tr = solution.extras.radiation_transmitted
    absorbed = solution.extras.volumetric_absorption
    isempty(tr) && return nothing
    I_max = maximum(tr .+ absorbed)
    I_max > 0 || return nothing
    frac = maximum(tr) / I_max
    if frac > 0.01
        @warn "Beer–Lambert radiation transmitted through the back face: up to " *
              "$(round(100 * frac, digits = 1)) % of the entering flux leaves the " *
              "domain (`solution.extras.radiation_transmitted`). If the back face " *
              "rests on an opaque holder, give the substrate an absorption coefficient."
    end
    return nothing
end

"""
	extract_solution_ale!(solution::PyrolysisSolution, raw_solution, problem, layout::StateLayout)

Extract results from ALE ODE solution into PyrolysisSolution structure.

For ALE solutions, node positions are part of the state vector, so we can
reconstruct the geometry at each saved time and compute volumes accordingly.
This ensures proper mass accounting with the evolving mesh.

With depletion handling the number of live cells at a saved time comes from
the recorded park events (`n_live_cells`); parked ghost slots of `T`, `ξ`
and `z` read as `NaN`, the thickness is the live thickness
`z[n_live+1] − z[1]`, and every sum runs over the live cells.

# Keyword Arguments
- `depletion::Union{Nothing,DepletionCallbackState}`: park-event record of
  the solve (`nothing` when depletion handling was off).
- `radiation_model::RadiationModel`: Radiation model for boundary flux diagnostics.
"""
function extract_solution_ale!(
    solution::PyrolysisSolution,
    raw_solution,
    problem,
    layout::StateLayout{NC};
    depletion::Union{Nothing,DepletionCallbackState} = nothing,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
) where {NC}
    mesh = problem.mesh
    material = problem.material
    NR = n_reactions(material)

    # Store raw solution for interpolation
    solution.raw_solution = raw_solution

    # Get return code
    solution.retcode = Symbol(raw_solution.retcode)

    times = raw_solution.t
    solution.t = copy(times)
    n_times = length(times)
    n_cells = layout.n_cells
    n_nodes = layout.n_nodes

    # Pre-allocate output arrays. Parked ghost slots read as NaN.
    solution.T = fill(NaN, n_cells, n_times)
    solution.ξ = fill(NaN, NC, n_cells, n_times)
    solution.z = fill(NaN, n_nodes, n_times)
    solution.surface_T = zeros(n_times)
    solution.mass_loss_rate = zeros(n_times)
    _allocate_extras!(solution, n_times)

    # Pre-allocate work arrays for reaction rate computation
    rates = zeros(NR)
    heats = zeros(NR)

    # Get boundary conditions for flux computation
    bc_set = problem.bc_set
    thermal_bc_top = get_thermal_bc(bc_set, Val(:top))
    thermal_bc_bottom = get_thermal_bc(bc_set, Val(:bottom))
    mass_bc_top = get_mass_bc(bc_set, Val(:top))
    mass_bc_bottom = get_mass_bc(bc_set, Val(:bottom))

    depletion === nothing || append!(solution.extras.park_times, depletion.park_times)

    # Extract data at each saved time
    for (out_idx, t) in enumerate(times)
        u = raw_solution.u[out_idx]
        n_live = depletion === nothing ? n_cells : n_live_cells(depletion, t)
        solution.extras.live_cells[out_idx] = n_live

        # Live temperatures (ghost slots stay NaN)
        T_profile = [get_temperature(u, layout, j) for j = 1:n_live]
        @inbounds for j = 1:n_live
            solution.T[j, out_idx] = T_profile[j]
        end

        # Live concentrations
        ξ_profile = zeros(NC, n_live)
        @inbounds for j = 1:n_live
            for k = 1:NC
                ξ_profile[k, j] = get_concentration(u, layout, j, k)
                solution.ξ[k, j, out_idx] = ξ_profile[k, j]
            end
        end

        # Live node positions and thickness
        z = get_node_positions(u, layout)
        @inbounds for j = 1:(n_live+1)
            solution.z[j, out_idx] = z[j]
        end
        thickness = z[n_live+1] - z[1]
        solution.extras.thickness_history[out_idx] = thickness

        # Compute cell volumes from current node positions and the state-derived A(t).
        A_t = cross_section_area_from_state(mesh, material, u, layout)
        volumes = compute_volumes_from_state(u, layout, A_t)

        # ================================================================
        # Solve for true boundary temperatures. The per-species gas mass flux
        # at each face is reconstructed the same way the live mass BC does,
        # so the blowing models see the flux the integrator saw.
        # ================================================================

        # Surface (top, boundary_id=1) — the live top cell
        T_cell_top = T_profile[n_live]
        ξ_top = ntuple(k -> ξ_profile[k, n_live], NC)
        k_top = effective_conductivity(ξ_top, T_cell_top, material)
        ε_top = effective_emissivity(ξ_top, T_cell_top, material)
        a_top = effective_absorptivity(ξ_top, T_cell_top, material)
        Δz_half_top = (z[n_live+1] - z[n_live]) / 2

        J_top = _surface_gas_flux(
            material,
            ξ_top,
            T_cell_top,
            mass_bc_top,
            Δz_half_top,
            t,
            has_darcy_flow(material),
        )

        state_top = solve_boundary_temperature(
            T_cell_top,
            k_top,
            Δz_half_top,
            thermal_bc_top,
            material,
            J_top,
            t,
            ε_top;
            a_eff = a_top,
            radiation_model = radiation_model,
        )
        solution.surface_T[out_idx] = state_top.T_surface
        solution.extras.surface_blowing_ratio[out_idx] = state_top.blowing_ratio
        solution.extras.surface_pooling_factor[out_idx] = state_top.pooling_factor
        solution.extras.surface_face_fraction[out_idx] = state_top.face_fraction
        solution.extras.surface_face_buoyancy[out_idx] = state_top.face_buoyancy

        # Back (bottom, boundary_id=2) — vents through the user-supplied bottom
        # mass BC (zero for the default impermeable holder; the RHS applies
        # a permeable back-face BC, so the MLR must count it).
        T_cell_bottom = T_profile[1]
        ξ_bottom = ntuple(k -> ξ_profile[k, 1], NC)
        k_bottom = effective_conductivity(ξ_bottom, T_cell_bottom, material)
        ε_bottom = effective_emissivity(ξ_bottom, T_cell_bottom, material)
        a_bottom = effective_absorptivity(ξ_bottom, T_cell_bottom, material)
        Δz_half_bottom = (z[2] - z[1]) / 2
        J_bottom = _surface_gas_flux(
            material,
            ξ_bottom,
            T_cell_bottom,
            mass_bc_bottom,
            Δz_half_bottom,
            t,
            has_darcy_flow(material),
            Val(:bottom),
        )
        state_bottom = solve_boundary_temperature(
            T_cell_bottom,
            k_bottom,
            Δz_half_bottom,
            thermal_bc_bottom,
            material,
            J_bottom,
            t,
            ε_bottom;
            a_eff = a_bottom,
            radiation_model = radiation_model,
        )
        solution.extras.back_temperature[out_idx] = state_bottom.T_surface

        # Total mass (integrate ξ × V over the live cells)
        total_m = 0.0
        for j = 1:n_live
            for k = 1:NC
                total_m += ξ_profile[k, j] * volumes[j]
            end
        end
        # Normalize by original area for cone-gauge-style mass loading [kg/m²]
        solution.extras.total_mass[out_idx] = total_m / mesh.cross_section_area_initial

        # ================================================================
        # Boundary flux computation (using solved surface temperature)
        # ================================================================

        # Surface cell properties - reuse variables from boundary temperature solve
        T_surface_solved = state_top.T_surface

        # Store effective properties
        solution.extras.surface_emissivity[out_idx] = ε_top
        solution.extras.surface_conductivity[out_idx] = k_top

        # Compute heat flux components using SOLVED surface temperature.
        # The face the Newton applied the blowing models to, so the convective
        # decomposition matches the corrected flux.
        heat_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_top,
            T_surface_solved,
            t,
            ε_top,
            radiation_model,
            BlowingFace(J_top, T_surface_solved, material);
            a_eff = a_top,
        )
        solution.extras.surface_heat_flux_absorbed[out_idx] = heat_fluxes.absorbed
        solution.extras.surface_heat_flux_reradiation[out_idx] = heat_fluxes.reradiation
        solution.extras.surface_heat_flux_convection[out_idx] = heat_fluxes.convection
        # Dirichlet: report the solved conducted flux.
        solution.extras.surface_heat_flux_net[out_idx] =
            thermal_bc_top isa DirichletThermalBC ? state_top.q_conducted : heat_fluxes.net

        # Conductive heat flux at surface: use solved boundary state
        # This is consistent with the boundary temperature solve and ensures energy balance
        solution.extras.surface_heat_flux_conduction[out_idx] = state_top.q_conducted

        # Volumetric radiative absorption (BEER_LAMBERT only) on the live ALE
        # geometry: Σ q'''·Δz = α_surf·q_ext·(1 − e^{−τ_tot}); see the
        # Eulerian extractor for the convention.
        if radiation_model == BEER_LAMBERT
            I_inc, has_explicit_α = extract_radiative_flux_from_bcs(bc_set, Val(:top), t)
            I_eff = has_explicit_α ? I_inc : ε_top * I_inc
            τ_tot = 0.0
            for j = 1:n_live
                ξ_j = ntuple(k -> ξ_profile[k, j], NC)
                τ_tot += effective_absorption(ξ_j, material) * (z[j+1] - z[j])
            end
            solution.extras.volumetric_absorption[out_idx] = I_eff * (1.0 - exp(-τ_tot))
            solution.extras.radiation_transmitted[out_idx] = I_eff * exp(-τ_tot)
        end

        # Surface MLR: per-species gas flux already used in the boundary
        # solve. `surface_mass_flux` is the top face only; the total MLR also
        # counts the bottom film/venting flux (zero for the default
        # impermeable holder).
        total_mass_flux = sum(J_top)
        bottom_mass_flux = sum(J_bottom)
        solution.extras.surface_mass_flux[out_idx] = total_mass_flux

        # Variable cross-section MLR diagnostics
        # total_mass_flux is the per-current-area surface flux [kg/(m²·s)] —
        # i.e., MLR_local. ṁ_total = MLR_local · A(t). MLR_gauge = ṁ_total / A_0.
        solution.extras.cross_section_area_history[out_idx] = A_t
        solution.mass_loss_rate[out_idx] = total_mass_flux + bottom_mass_flux
        solution.extras.MLR_gauge[out_idx] =
            (total_mass_flux + bottom_mass_flux) * A_t / mesh.cross_section_area_initial

        # Back face heat flux (boundary_id=2, z=0) - use solved boundary temperature
        T_back_solved = state_bottom.T_surface
        back_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_bottom,
            T_back_solved,
            t,
            ε_bottom,
            radiation_model,
            BlowingFace(J_bottom, T_back_solved, material);
            a_eff = a_bottom,
        )
        solution.extras.back_heat_flux[out_idx] =
            thermal_bc_bottom isa DirichletThermalBC ? state_bottom.q_conducted :
            back_fluxes.net

        # Gas generation rate: compute from reaction rates (physics-based)
        # This is NOT the actual mass loss rate (surface flux), but rather the
        # volumetric production of gases throughout the sample.
        gen_rate = 0.0
        for j = 1:n_live
            T = T_profile[j]
            ξ = ntuple(k -> ξ_profile[k, j], NC)
            cell_vol = volumes[j]

            # Evaluate reaction rates at this cell
            Physics.evaluate_reactions!(rates, heats, ξ, T, material)

            # Sum gas production from each reaction
            for r_idx = 1:NR
                rxn = material.reactions[r_idx]
                r = rates[r_idx]

                # Gas production: sum mass yields for any gas-phase product
                for (p_idx, ν) in rxn.products
                    if material.components[p_idx] isa GaseousComponent
                        gen_rate += ν * r * cell_vol
                    end
                end
            end
        end
        # Normalize by original area for cone-gauge-style generation rate [kg/(m²·s)]
        solution.extras.gas_generation_rate[out_idx] =
            gen_rate / mesh.cross_section_area_initial
    end

    _warn_transmitted_radiation(solution, radiation_model)
    return solution
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
	get_profile(solution::PyrolysisSolution, field::Symbol, time_idx::Int)

Get a field profile at a specific saved time index.

# Arguments
- `solution`: The solution
- `field`: :T for temperature, :ξ for all concentrations, or :ξ1, :ξ2, etc.
- `time_idx`: Index into solution.t

# Returns
Vector of values at cell centers.
"""
function get_profile(solution::PyrolysisSolution, field::Symbol, time_idx::Int)
    if field == :T
        return solution.T[:, time_idx]
    elseif field == :ξ
        return solution.ξ[:, :, time_idx]
    else
        # Parse :ξ1, :ξ2, etc.
        field_str = string(field)
        if startswith(field_str, "ξ") || startswith(field_str, "xi")
            idx_str = replace(field_str, r"^(ξ|xi)" => "")
            comp_idx = parse(Int, idx_str)
            return solution.ξ[comp_idx, :, time_idx]
        end
    end
    error("Unknown field: $field. Use :T, :ξ, :ξ1, :ξ2, etc.")
end

"""
	get_cell_positions(solution::PyrolysisSolution)

Get cell center positions from the solution's mesh.
"""
function get_cell_positions(solution::PyrolysisSolution)
    return solution.problem.mesh.cell_centers
end

"""
	interpolate_solution(solution::PyrolysisSolution, t::Real)

Interpolate the full state vector at time t.
"""
function interpolate_solution(solution::PyrolysisSolution, t::Real)
    if solution.raw_solution === nothing
        error("Raw solution not available for interpolation")
    end
    return solution.raw_solution(t)
end
