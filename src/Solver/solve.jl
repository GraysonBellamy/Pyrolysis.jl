"""
Main solve interface for PyrolysisProblem.
"""

using ..Problem:
    PyrolysisProblem,
    PyrolysisSolution,
    validate_problem,
    apply_initial_conditions!,
    to_problem_def,
    build_workspace,
    StateLayout,
    state_length,
    create_state_vector,
    get_temperature,
    get_concentration,
    get_node_position,
    get_node_positions,
    compute_volumes_from_state,
    get_thickness_from_state,
    z_offset,
    chi_range,
    _allocate_extras!
using ..Residual: residual!, cross_section_area_from_state
using ..Geometry: cell_center
using ..Physics: evaluate_reactions!, ideal_gas_pressure
using ..Properties:
    RadiationCache,
    effective_emissivity,
    effective_absorption,
    effective_conductivity,
    effective_gas_transfer,
    effective_permeability,
    effective_gas_viscosity,
    porosity
using ..Discretization:
    compute_boundary_heat_fluxes,
    compute_surface_mass_flux,
    solve_boundary_temperature,
    boundary_darcy_velocity
using ..Materials: AbstractThermalBC, AbstractRadiationBC, GAS, T_REF, has_darcy_flow
using ..BoundaryConditions:
    get_radiative_flux_bc,
    extract_radiative_flux_from_bcs,
    has_radiative_flux_bc,
    RadiativeBC,
    CombinedThermalBC,
    DirichletThermalBC,
    get_thermal_bc,
    get_mass_bc,
    ambient_temperature
using LinearSolve: LUFactorization, KLUFactorization, SparspakFactorization  # LU for MKL workaround; Sparspak is the default sparse direct solver.
using SparseArrays: nnz  # For logging sparse Jacobian info

import CommonSolve: solve  # Extend CommonSolve.solve for our problem type

# ============================================================================
# Post-processing helpers (private)
# ============================================================================
#
# The saved-output extraction needs the same per-species gas mass flux at the
# boundary that the live RHS computes, so the boundary-temperature solve here
# matches what the integrator saw during the run. Without this, the saved
# `surface_temperature` and `surface_heat_flux_*` fields silently use a
# zero-blowing approximation while the live solve uses the full transpiration
# BC, producing diagnostic values that are inconsistent with the integrator
# state. These helpers also let us audit the transpiration energy sink as a
# stored field on `PyrolysisSolution`.
# (Long-term these belong in Discretization; for now they live alongside
# the extraction code that uses them.)

# Per-species outward gas mass flux at a boundary face [kg/(m²·s)],
# evaluated the same way the live mass BC application does: diffusive film
# law plus, in Darcy-Fick mode, the pressure-BC advective venting
# v_out·ξ/φ (same `boundary_darcy_velocity` kernel as the RHS; outflow only —
# inflow advects the zero-concentration ambient). `tag` is the boundary
# orientation: at the bottom the outward normal is −z, so outflow is v_z < 0
# (the bottom face vents too when a permeable mass BC is configured there;
# "positive = leaving the slab" at both boundaries).
@inline function _surface_gas_flux(
    material::Material{NC},
    ξ_face::NTuple{NC,Float64},
    T_face::Float64,
    mass_bc,
    d_half::Float64,
    t::Float64,
    darcy_active::Bool,
    ::Val{tag} = Val(:top),
) where {NC,tag}
    ρ = sum(ξ_face)
    λρ = effective_gas_transfer(ξ_face, T_face, material) * ρ
    # Single assignment (no conditional reassignment) so the `ntuple` closure
    # below captures a non-boxed value — avoids a JET possible-undef report.
    v_out_over_φ = if darcy_active
        φ = porosity(ξ_face, T_face, material)
        P_int = ideal_gas_pressure(ξ_face, T_face, φ, material)
        κ_int = effective_permeability(ξ_face, T_face, material)
        μ_int = effective_gas_viscosity(T_face)
        v_z = boundary_darcy_velocity(mass_bc, P_int, κ_int, μ_int, d_half, Val(tag))
        v_out = tag === :top ? v_z : -v_z
        v_out > 0.0 ? v_out / max(φ, 1e-12) : 0.0
    else
        0.0
    end
    ntuple(Val(NC)) do k
        comp = material.components[k]
        if is_gas(comp) && ξ_face[k] > 0.0
            compute_surface_mass_flux(mass_bc, ξ_face[k], λρ, ρ, d_half, t, k) +
            v_out_over_φ * ξ_face[k]
        else
            0.0
        end
    end
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
- P1_QUASI_STEADY is requested on the core solver path (use Experimental)
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

    # P1 lives outside the core unified residual.
    if model == P1_QUASI_STEADY
        @warn """
          P1_QUASI_STEADY is experimental and is not wired into the core unified
          solve path. Use Pyrolysis.Experimental.RadiationP1 primitives directly.
          """ maxlog = 1
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
	_validate_solve_kwargs(; extra_kwargs, jacobian, handle_depletion,
	                          radiation_model, energy_form, use_transpiration_bc,
	                          min_thickness, depletion_threshold, min_cells,
	                          diagnostics_stride, verify_jacobian)

Catch every cross-kwarg compatibility error at solve() entry, before any
allocation. Individual numeric or symbol checks belong here too whenever they
can be evaluated from kwargs alone (no context).
"""
function _validate_solve_kwargs(;
    extra_kwargs,
    jacobian,
    handle_depletion,
    radiation_model,
    energy_form,
    use_transpiration_bc,
    min_thickness,
    depletion_threshold,
    min_cells,
    diagnostics_stride,
    verify_jacobian,
)
    # Unknown kwarg rejection.
    _reject_unknown_kwargs(extra_kwargs)

    # Cross-kwarg combinations.
    if jacobian !== nothing && handle_depletion
        error(
            "solve(): the `jacobian` kwarg is incompatible with " *
            "`handle_depletion = true` (state-vector resizing is not " *
            "supported by the pluggable backends). Drop one of the two.",
        )
    end
    if radiation_model == P1_QUASI_STEADY
        error(
            "solve(): P1_QUASI_STEADY radiation lives in " *
            "`Pyrolysis.Experimental.RadiationP1` and is not part of the " *
            "core unified residual.",
        )
    end

    # Single-kwarg shape checks that fail fast.
    if !(energy_form === :standard || energy_form === :with_generation_sink)
        throw(
            ArgumentError(
                "solve(): energy_form must be :standard or " *
                ":with_generation_sink, got $(repr(energy_form)).",
            ),
        )
    end
    # Guard the gas-convective-enthalpy double-count (see to_problem_def): the
    # volumetric S_conv is always active, so enabling the surface transpiration
    # term too over-applies the gas enthalpy sink by ∫S_conv.
    if use_transpiration_bc
        throw(
            ArgumentError(
                "solve(): use_transpiration_bc=true double-counts gas convective " *
                "enthalpy (S_conv is always active and overlaps the surface " *
                "transpiration term by ∫S_conv). Use use_transpiration_bc=false " *
                "(the default, S_conv-only ThermaKin Eq.17 convention).",
            ),
        )
    end
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
# When `jacobian = nothing`, pick `Structured(strategy = DirectSolve())`
# for both Eulerian and ALE problems. ALE callers wanting the
# Schur-style approximation can pass
# `JacobianSpec(Structured(strategy = ApproxSolve()), problem)`.
#
# ApproxSolve is deliberately NOT the ALE default despite being ~8-10×
# faster end-to-end on smooth horizons with identical integrator stats
# (perf fixtures, 2026-06-11: t = 5 s and 15 s gave the same naccept/njacs,
# nf within 4%). Near burnout the dropped ALE couplings cost real
# robustness: on the ale_beer_lambert_darcy fixture at t = 30 s (violent
# depletion tail; 1681 accepted steps) ApproxSolve's larger Newton error
# admits far-from-physical trial states — it formerly aborted with a
# DomainError from a negative-temperature iterate reaching Sutherland's
# T^1.5. The gas-property kernels now clamp (T_GAS_PROPERTY_FLOOR,
# src/Properties/effective_properties.jl) and that solve completes, but
# only by surviving wilder iterates than DirectSolve ever produces.
# Use ApproxSolve explicitly when the simulated horizon stays away from
# full surface depletion.

function _default_jacobian_spec(def, ws, use_ale::Bool, handle_depletion::Bool)
    if handle_depletion
        # Depletion resizes the state vector mid-solve via `resize!`; the
        # pluggable backends carry an immutable `JacobianSpec` cache built
        # at problem-setup time that does not survive resize. No Jacobian
        # is supplied — the integrator falls back to its own colored FD.
        return nothing
    end
    return JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
end

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
    sample_times = collect(get(verify_kwargs, :sample_times, [0.0]))
    # Build the dense reference once and reuse it across firings.
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)

    condition =
        (u, t, integrator) -> begin
            # Fire once we cross any sample time. Tolerance: half the
            # integrator's current dt.
            dt = max(integrator.dt, eps(Float64))
            any(abs(t - st) <= 0.5 * dt for st in sample_times)
        end

    affect! =
        (integrator) -> begin
            report = verify_jacobian_at(
                spec_active,
                def,
                ws,
                integrator.u,
                integrator.t;
                reference_spec = ref_spec,
                rtol = rtol,
                atol = atol,
            )
            if !report.passed
                @warn "verify_jacobian: backend $(typeof(spec_active.backend)) " *
                      "deviates from DenseAD reference at t=$(integrator.t): " *
                      "max_abs_dev=$(report.max_abs_dev), max_rel_dev=$(report.max_rel_dev), " *
                      "worst_index=$(report.worst_index)"
            end
            return nothing
        end

    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# ============================================================================
# Linear solver / integrator selection
# ============================================================================

_is_sparse_linsolve(ls) = ls isa SparspakFactorization || ls isa KLUFactorization

"""
	_resolve_integrator(; integrator, linear_solver, needs_dense_linsolve=false)

Resolve the user's integrator / linear-solver request against the Jacobian
path. `needs_dense_linsolve = true` means no sparse `jac_prototype` will be
attached to the `ODEFunction` (the `jacobian === nothing` colored-FD path —
e.g. under `handle_depletion` — or a dense Jacobian backend), where a sparse
factorization cannot run.

The user's integrator choice is always preserved: an incompatible sparse
factorization is swapped for `LUFactorization` **within** the user's
algorithm (constructor re-call, or a rebuild of a fully constructed
algorithm keeping its type and `autodiff`), with a warning — never by
silently substituting a different integrator.
"""
function _resolve_integrator(;
    integrator,
    linear_solver,
    needs_dense_linsolve::Bool = false,
)
    default_ls = needs_dense_linsolve ? LUFactorization() : SparspakFactorization()
    ls = linear_solver === nothing ? default_ls : linear_solver
    if needs_dense_linsolve && _is_sparse_linsolve(ls)
        @warn "linear_solver = $(typeof(ls)) requires a sparse Jacobian " *
              "prototype, but this solve has none (jacobian = nothing / " *
              "depletion colored-FD path, or a dense Jacobian backend). " *
              "Using LUFactorization instead."
        ls = LUFactorization()
    end

    if integrator === nothing
        return KenCarp4(autodiff = AutoFiniteDiff(), linsolve = ls)
    elseif integrator isa Type || integrator isa Function
        return integrator(autodiff = AutoFiniteDiff(), linsolve = ls)
    elseif linear_solver === nothing
        # Fully constructed algorithm.
        alg = integrator
        if needs_dense_linsolve &&
           hasproperty(alg, :linsolve) &&
           _is_sparse_linsolve(alg.linsolve)
            alg_ctor = Base.typename(typeof(alg)).wrapper
            @warn "integrator = $(alg_ctor)(linsolve = $(typeof(alg.linsolve))) " *
                  "requires a sparse Jacobian prototype, but this solve has " *
                  "none. Rebuilding $(alg_ctor) with LUFactorization — the " *
                  "algorithm type and autodiff are preserved; any other " *
                  "custom algorithm options are reset to their defaults."
            return if hasproperty(alg, :autodiff)
                alg_ctor(autodiff = alg.autodiff, linsolve = LUFactorization())
            else
                alg_ctor(linsolve = LUFactorization())
            end
        end
        return alg
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
- `problem`: The PyrolysisProblem to solve

# Keyword Arguments
- `integrator=nothing`: ODE integrator constructor (e.g. `KenCarp4`) or fully
  constructed algorithm (e.g. `KenCarp4(linsolve=SparspakFactorization())`).
  `nothing` chooses `KenCarp4` with the default `linear_solver`.
- `linear_solver=nothing`: LinearSolve backend. `nothing` chooses
  `SparspakFactorization()` for sparse Jacobians. Pair with `integrator` only
  when `integrator` is a constructor.
- `jacobian=nothing`: Jacobian handle (`JacobianSpec`). `nothing`
  picks `Structured(strategy = DirectSolve())`. Pass an explicit
  `JacobianSpec(Structured(strategy = ApproxSolve()), problem)` for the
  Schur-style approximation that drops the structured ALE couplings.
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
- `min_thickness=1e-6`: Minimum thickness before termination (ALE only) [m]
- `handle_depletion=false`: Enable surface cell depletion handling (ALE only)
- `depletion_threshold=0.05`: Thickness fraction threshold for cell merging
  (0.05 = 5% of the cell's own initial thickness)
- `min_cells=2`: Minimum number of cells before terminating instead of merging
- `use_transpiration_bc=false`: Opt-in surface gas sensible-enthalpy transpiration. **Off by
  default** — the volumetric advection source S_conv already carries the gas enthalpy
  (ThermaKin Eq.17 convention); enabling this alongside S_conv double-counts and is rejected
  at problem construction. Retained only for closed-ledger experiments with S_conv disabled.
- `use_blowing_correction=false`: Opt-in Couette-film blowing correction to every
  convective BC in the surface energy balance: `h → h·B/(e^B − 1)` with
  `B = ṁ″·c_p,g/h` built from the outward gas mass flux at the face. The
  momentum-side companion of the transpiration term — unlike it, fully compatible
  with the S_conv-only convention (no enthalpy double-count). Suppresses convective
  cooling at peak burning (B = 4 → h_eff/h = 0.075); the back face is unaffected
  automatically (no outflow → B = 0). Calibrated examples should enable this only
  together with flame heat-flux feedback (see the physics roadmap).
- `warm_start=nothing`: Reuse a prior `PyrolysisSolution`'s final state as the initial state
  (`u0`). Component count, mesh size, and ALE/χ modes must match. `tspan[1]` still defines
  the integrator's start time, so the warm-started solve treats the prior final state as
  fresh initial data.

# Jacobian backend

Pass `jacobian = JacobianSpec(backend, problem)` to override the
auto-default. The production backend is `Structured`; behavior is
selected by its strategy:

| Spec                                                   | When to use                              |
|--------------------------------------------------------|------------------------------------------|
| `Structured(strategy = DirectSolve())`                 | Default; exact assembly + sparse direct  |
| `Structured(strategy = StructuredSolve())`             | Apply structured couplings without densify |
| `Structured(strategy = ApproxSolve())`                 | Drops structured couplings (ALE Schur-style approximation) |
| `DenseAD_ForwardDiff()` (verifier reference only)      | Internal use — reference for `verify_jacobian` |

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
When `handle_depletion=true` (requires `use_ale=true`), the solver monitors
every cell's thickness. When a cell shrinks below `depletion_threshold`
fraction of its **own** initial thickness (a per-cell reference, so stretched
meshes are judged cell-by-cell), it is merged with a neighbor. This prevents
the ρcₚ→0 singularity that causes temperature blow-up in materials that fully
decompose. The mesh is dynamically resized during integration.

# Recommended Solvers
The default is `KenCarp4(autodiff = AutoFiniteDiff(),
linsolve = SparspakFactorization())` with
`Structured(strategy = DirectSolve())`. This is the production baseline for
stiff pyrolysis runs.

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
With `jacobian = nothing` (the default), every mode uses
`Structured(strategy = DirectSolve())`. The `Structured` backend
assembles per-operator contributions into a `StructuredJacobian` and
densifies any structured couplings (ALE prefix-sum, χ rank-1, Beer-
Lambert suffix-sum) at solve time so the linear-solve sees a true
sparse Jacobian. For ALE problems where the dense cumulative C block
dominates factorization cost, opt into the Schur-style approximation
via `Structured(strategy = ApproxSolve())` — this drops the structured
couplings, at the cost of an approximate Newton search direction.

# Linear Solver
Sparspak is the default sparse direct factorization for pyrolysis problems.
Pass `linear_solver = KLUFactorization()`, another LinearSolve backend, or a
fully constructed `integrator` when you need to compare factorization choices.

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

# With ALE — auto-default Jacobian is Structured(strategy=DirectSolve()):
solution = solve(problem, use_ale=true,
    integrator=KenCarp4(linsolve=SparspakFactorization()))

# Opt into the ALE Schur-style approximation:
solution = solve(problem, use_ale=true,
    jacobian=JacobianSpec(Structured(strategy=ApproxSolve()), problem))

# With surface cell depletion handling (prevents blow-up):
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
    # `jacobian` is the single Jacobian handle. `nothing` →
    # `Structured(strategy = DirectSolve())`.
    jacobian::Union{JacobianSpec,Nothing} = nothing,
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
    # First-class runtime conservation diagnostics.
    diagnostics::Bool = false,
    diagnostics_stride::Int = 10,
    # `energy_form` selects the bulk energy formulation:
    #   :standard               — Formulation A (drops the gas-generation sink)
    #   :with_generation_sink  — adds the −Σ Δh_g·∇·J_g cell source. Requires
    #                             every reaction heat to be referenced to the
    #                             T_ref (≈298.15 K) datum; with DSC-style heats
    #                             quoted at the reaction temperature it
    #                             double-counts the gas sensible enthalpy
    #                             (see `Problem.to_problem_def`, which warns).
    # The :standard gas-energy convention (S_conv only) matches ThermaKin Eq.17.
    energy_form::Symbol = :standard,
    # Opt-in surface gas transpiration enthalpy. OFF by default: the volumetric
    # S_conv already carries gas enthalpy, and combining the two double-counts
    # (rejected at problem construction). Mass transport remains active either way.
    use_transpiration_bc::Bool = false,
    # Opt-in Couette-film blowing correction h → h·B/(e^B − 1) on convective BCs
    # in the surface energy balance (B = ṁ″·c_p,g/h from the outward gas flux).
    use_blowing_correction::Bool = false,

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
        jacobian = jacobian,
        handle_depletion = handle_depletion,
        radiation_model = radiation_model,
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        min_thickness = min_thickness,
        depletion_threshold = depletion_threshold,
        min_cells = min_cells,
        diagnostics_stride = diagnostics_stride,
        verify_jacobian = verify_jacobian,
    )

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

    # Apply initial conditions to mesh
    apply_initial_conditions!(problem)

    # Create solution container
    solution = PyrolysisSolution(problem)

    # Start timing
    start_time = time()

    # Get radiation flux function from boundary conditions
    radiation_flux = _get_radiation_flux_function(problem)

    # Decide how the external radiation flux feeds the ODE.
    # - P1_QUASI_STEADY: handled via Marshak BC, no I_radiation here.
    # - NO_RADIATION:    no flux at all.
    # - SURFACE_ABSORPTION / BEER_LAMBERT: the flux is re-sampled at EVERY RHS
    #   evaluation (the old two-point t=0/t=1 equality
    #   detection silently froze any flux that merely *happened* to agree at
    #   those instants — e.g. q₀(1+sin 2πt) or a post-1 s ramp — at I(0) for
    #   the volumetric deposit). One user-function call per RHS is negligible
    #   against the residual; constant fluxes lose nothing.
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

    def = to_problem_def(
        problem;
        use_ale = use_ale,
        radiation_model = radiation_model,
        radiation_intensity = I_radiation,
        radiation_kirchhoff = radiation_kirchhoff,
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        use_blowing_correction = use_blowing_correction,
    )
    ws = build_workspace(problem, def; min_thickness = min_thickness)
    u0 = create_state_vector(ws.layout, ws.mesh)
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

    # Build the ODE problem from the Jacobian spec. The user either supplies
    # a `JacobianSpec` via the `jacobian` kwarg, or accepts the auto-default
    # (`Structured(strategy = DirectSolve())`).
    local ode_problem
    local jac_spec::Union{Nothing,JacobianSpec} =
        jacobian === nothing ? _default_jacobian_spec(def, ws, use_ale, handle_depletion) :
        jacobian

    # The Jacobian path fixes the linear-solver family: a sparse prototype
    # wants a sparse factorization, while the no-spec colored-FD path and the
    # dense backends need a dense one. Resolve the integrator only now so the
    # default linsolve matches the path, and an incompatible explicit choice
    # is rebuilt inside the *user's* algorithm (with a warning) instead of
    # being silently replaced by a hard-coded KenCarp4.
    needs_dense_linsolve = jac_spec === nothing || jac_spec.sparsity_phys === nothing
    solver_alg = _resolve_integrator(; integrator, linear_solver, needs_dense_linsolve)

    if jac_spec === nothing
        # No Jacobian: fall back to the integrator's own colored FD. Used
        # under depletion or by users opting out explicitly via
        # `jacobian = nothing`.
        ode_problem = ODEProblem(rhs!, u0, tspan)
    else
        # Build the ODEFunction from the active spec. Backends with
        # sparsity provide a `jac_prototype` for sparse LU factorization;
        # the verifier-only `DenseAD_ForwardDiff` does not.
        jac_fn! = let s = jac_spec, def = def, ws = ws
            (J, u, p, t) -> compute_jacobian!(J, s, u, def, t, ws)
        end

        jac_prototype = if jac_spec.sparsity_phys !== nothing
            # Build the union of physics + geometry sparsity as the
            # allocation template. Both are `SparseMatrixCSC{Bool,Int}`;
            # the Float64 working matrix is materialized by OrdinaryDiffEq.
            S_geom = jac_spec.sparsity_geom
            S_union =
                S_geom === nothing ? jac_spec.sparsity_phys :
                (jac_spec.sparsity_phys .| S_geom)
            # Convert Bool sparsity to Float64 prototype.
            Float64.(S_union)
        else
            nothing
        end

        ode_fn_with_jac = if jac_prototype !== nothing
            OrdinaryDiffEq.ODEFunction(rhs!; jac = jac_fn!, jac_prototype = jac_prototype)
        else
            OrdinaryDiffEq.ODEFunction(rhs!; jac = jac_fn!)
        end
        ode_problem = ODEProblem(ode_fn_with_jac, u0, tspan)

        # Log + warn about linsolve choice. Sparse Jacobians want a
        # sparse linsolve (Sparspak is 28% faster than KLU for the
        # block-tridiagonal pattern); dense Jacobians require a dense one.
        if jac_prototype !== nothing
            @info "Using `$(typeof(jac_spec.backend))` Jacobian " *
                  "with sparse storage ($(nnz(jac_prototype)) non-zeros)"
            if hasproperty(solver_alg, :linsolve) &&
               solver_alg.linsolve isa ILUGMRESFactorization
                @warn "Using ILU-GMRES iterative solver (τ=$(solver_alg.linsolve.τ)). " *
                      "Note: Sparspak is 1.7-3.3× faster for production pyrolysis " *
                      "simulations due to Newton step rejection from inexact GMRES solves. " *
                      "Consider KenCarp4(linsolve=SparspakFactorization()) instead." maxlog =
                    1
            elseif hasproperty(solver_alg, :linsolve) &&
                   !(solver_alg.linsolve isa SparspakFactorization) &&
                   !(solver_alg.linsolve isa KLUFactorization)
                @warn "Sparse Jacobian — Sparspak is 5× faster than the default " *
                      "dense solver. Pass linear_solver=SparspakFactorization() " *
                      "for best performance."
            end
        end
        # Dense backends (no jac_prototype) already received a dense linsolve
        # in _resolve_integrator above.
    end

    # ============================================================================
    # MKL LU Factorization Workaround for Depletion Handling
    # ============================================================================
    # When using handle_depletion=true, the state vector is dynamically resized
    # during integration via resize!(integrator, new_length). For large problems
    # (state vector length > ~500), LinearSolve.jl's DefaultLinearSolver switches
    # to MKLLUFactorization which has a bug where the pivot array (ipiv) is not
    # properly resized, causing "DimensionMismatch: ipiv has length X, but needs
    # to be Y" errors.
    #
    # Workaround: When depletion handling is enabled and the state vector is large,
    # we switch to Julia's native LUFactorization which handles resize! correctly.
    # KLU also handles resize! correctly, so this only applies when using the
    # default linear solver (no explicit linsolve set).
    #
    # Reference: https://github.com/SciML/LinearSolve.jl (potential upstream bug)
    # ============================================================================
    # ============================================================================
    # MKL LU Factorization Workaround for Depletion Handling
    # ============================================================================
    # When using handle_depletion=true with resize!, MKL's pivot array may not
    # resize correctly. KLU and LUFactorization handle resize! properly.
    # Since we now default to KLU, this only warns if using the default linsolve.
    MKL_RESIZE_THRESHOLD = 500
    if use_ale && handle_depletion && length(u0) > MKL_RESIZE_THRESHOLD
        if hasproperty(solver_alg, :linsolve) && solver_alg.linsolve === nothing
            @warn "Large state vector ($(length(u0)) elements) with depletion handling. " *
                  "The default MKL linear solver may fail on resize!. " *
                  "Pass linear_solver=SparspakFactorization() to avoid this."
        end
    end

    # Build solver options
    solve_kwargs = Dict{Symbol,Any}(
        :abstol => abstol,
        :reltol => reltol,
        :maxiters => maxiters,
        :dtmin => dtmin,
        :dtmax => dtmax,
        :dt => problem.dt_initial,
        :progress => false,  # Disable OrdinaryDiffEq's progress, we use our own
    )
    if sensealg !== nothing
        solve_kwargs[:sensealg] = sensealg
    end

    # When handling depletion with resize!, we cannot use saveat with interpolation
    # because the interpolation doesn't handle state vector size changes.
    # Instead, we save every step and filter to desired times in post-processing.
    if use_ale && handle_depletion
        # Save every step (dense output without interpolation)
        solve_kwargs[:save_everystep] = true
        # Store the desired saveat times for post-processing
        requested_saveat = saveat
    else
        solve_kwargs[:saveat] = saveat
        requested_saveat = nothing
    end

    # Build callback list
    callbacks_list = []

    # StepUpdateCallback — mirror integrator state into the mutable mesh at the
    # end of each accepted step. Residual / Jacobian assembly
    # is mesh-read-only, so we need an explicit sync point for AMR, depletion
    # checks, and post-processing readers. Push first so it runs before any
    # downstream callback that reads mesh.
    push!(callbacks_list, create_step_update_callback(ws, material))

    # Create shared callback state if using depletion handling (so both callbacks see updates)
    depletion_cb_state = nothing
    if use_ale && handle_depletion
        # Pass mesh so DepletionCallbackState can compute per-cell initial thicknesses
        depletion_cb_state = DepletionCallbackState(ws.layout, mesh)
    end

    # Add thickness termination callback for ALE mode
    # Use shared state if available so it sees updated sv after cell merges
    if use_ale
        sv_ref = depletion_cb_state !== nothing ? depletion_cb_state : ws.layout
        thickness_cb = create_thickness_termination_callback(sv_ref, Float64(min_thickness))
        push!(callbacks_list, thickness_cb)
    end

    # Add surface cell depletion callback for ALE mode
    if use_ale && handle_depletion
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

    # Add progress callback if requested
    prog = nothing
    if progress
        t_start, t_end = tspan
        t_span_length = t_end - t_start
        last_progress = Ref(0)

        # Create progress bar
        prog = Progress(100; desc = "Solving: ", showspeed = true, enabled = true)

        # Use PeriodicCallback to update progress bar at regular time intervals
        check_interval = t_span_length / 100  # Check every 1% of simulation time

        progress_cb = PeriodicCallback(
            integrator -> begin
                t = integrator.t
                current_progress = floor(Int, 100 * (t - t_start) / t_span_length)
                if current_progress > last_progress[]
                    update!(prog, current_progress)
                    last_progress[] = current_progress
                end
            end,
            check_interval;
            initial_affect = false,
            save_positions = (false, false),  # Don't save solution at progress check times
        )
        push!(callbacks_list, progress_cb)
    end

    # Add diagnostics callback if requested. Snapshot mode — records
    # E_matrix, E_total, gas_storage_gap at stride boundaries. Cumulative
    # ledger fields (discrete/physical residuals) are filled by the
    # unified residual layer.
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

                integrator -> begin
                    stride_count[] += 1
                    if stride_count[] % stride == 0
                        # Sync mesh state from the integrator's u so the trackers
                        # read the current cell states.
                        update_mesh_from_state!(mesh, integrator.u, layout)
                        push!(
                            log,
                            snapshot_diagnostics_row(Float64(integrator.t), mesh, material),
                        )
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
    if verify_jacobian !== nothing && jac_spec !== nothing
        push!(callbacks_list, _build_verify_callback(jac_spec, def, ws, verify_jacobian))
    elseif verify_jacobian !== nothing
        @warn "verify_jacobian was requested but no Jacobian spec is " *
              "active (handle_depletion=true falls back to integrator FD); " *
              "verification skipped."
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

    # Finish progress bar if it was created
    if prog !== nothing
        finish!(prog)
    end

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
            requested_times = requested_saveat,
            radiation_model = radiation_model,
            use_transpiration_bc = use_transpiration_bc,
            use_blowing_correction = use_blowing_correction,
        )
    else
        extract_solution!(
            solution,
            raw_solution,
            problem;
            radiation_model = radiation_model,
            use_transpiration_bc = use_transpiration_bc,
            use_blowing_correction = use_blowing_correction,
        )
    end

    return solution
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
	_apply_warm_start!(u0, problem, layout, prev::PyrolysisSolution)

Overwrite `u0` with the final raw state from a prior `PyrolysisSolution`.

Compatibility checks:
- `prev.raw_solution` must exist (the prior solve must have populated it).
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

Note: The radiation_model passed to solve() determines how the flux is used:
- SURFACE_ABSORPTION: flux enters via the surface energy balance (boundary Newton)
- BEER_LAMBERT: flux used for volumetric Beer-Lambert absorption (I_radiation)
- P1_QUASI_STEADY: experimental; not routed through core solve
"""
function _get_radiation_flux_function(problem::PyrolysisProblem)
    if has_radiative_flux_bc(problem.bc_set, Val(:top))
        # Return a function that computes absorbed flux
        # This is used by BEER_LAMBERT and P1 models for volumetric absorption
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
	_extract_radiative_params(bc::AbstractThermalBC) -> Tuple{Any, Any}

Extract emissivity and T_ambient from a (possibly combined) thermal BC.
Returns (emissivity, T_ambient) with defaults (0.9, 300.0) if not found.
"""
_extract_radiative_params(bc::AbstractThermalBC) = (0.9, 300.0)

function _extract_radiative_params(bc::RadiativeBC)
    return (bc.ε, bc.T_inf)
end

function _extract_radiative_params(bc::CombinedThermalBC)
    for inner_bc in bc.bcs
        if inner_bc isa RadiativeBC
            return (inner_bc.ε, inner_bc.T_inf)
        end
    end
    # Recursively check for nested combined BCs
    for inner_bc in bc.bcs
        if inner_bc isa CombinedThermalBC
            result = _extract_radiative_params(inner_bc)
            if result[1] != 0.9 || result[2] != 300.0
                return result
            end
        end
    end
    return (0.9, 300.0)
end

"""
	extract_solution!(solution::PyrolysisSolution, raw_solution, problem; radiation_model)

Extract results from raw ODE solution into PyrolysisSolution structure.
"""
function extract_solution!(
    solution::PyrolysisSolution,
    raw_solution,
    problem;
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    use_transpiration_bc::Bool = false,
    use_blowing_correction::Bool = false,
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

    # Pre-allocate output arrays
    n_times = length(times)
    solution.T = Matrix{Float64}(undef, n_cells_val, n_times)
    solution.ξ = Array{Float64,3}(undef, NC, n_cells_val, n_times)
    solution.z = Matrix{Float64}(undef, 0, 0)             # Eulerian: no node positions
    solution.surface_T = zeros(n_times)
    solution.mass_loss_rate = zeros(n_times)
    solution.energy_residual = zeros(n_times)
    _allocate_extras!(solution, n_times)

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

        # Extract temperatures (write directly into T[:, i])
        T_profile = zeros(n_cells_val)
        @inbounds for j = 1:n_cells_val
            T_profile[j] = u[j]
            solution.T[j, i] = T_profile[j]
        end

        # Extract concentrations (write directly into ξ[:, :, i])
        ξ_profile = zeros(NC, n_cells_val)
        @inbounds for j = 1:n_cells_val
            for k = 1:NC
                ξ_profile[k, j] = u[n_cells_val+(k-1)*n_cells_val+j]
                solution.ξ[k, j, i] = ξ_profile[k, j]
            end
        end

        # ================================================================
        # Solve for true boundary temperatures (not cell-center values).
        # We reconstruct the per-species gas mass flux at the top face the
        # same way the live mass BC application does, then pass it to the
        # boundary-temperature solve so transpiration cooling is included
        # — keeping the saved diagnostics consistent with the integrator
        # state. The back face has no mass BC (insulated sample holder),
        # so J_bottom is zero by construction.
        # ================================================================
        T_ref_top = ambient_temperature(thermal_bc_top, t, Float64(T_REF))
        T_ref_bottom = ambient_temperature(thermal_bc_bottom, t, Float64(T_REF))

        # Surface (top, boundary_id=1) — geometry, properties, gas flux
        T_cell_top = T_profile[n_cells_val]
        ξ_top = ntuple(k -> ξ_profile[k, n_cells_val], NC)
        k_top = effective_conductivity(ξ_top, T_cell_top, material)
        ε_top = effective_emissivity(ξ_top, T_cell_top, material)
        Δz_half_top =
            mesh.cell_volumes[n_cells_val] /
            mesh.face_areas[mesh.cells[n_cells_val].face_ids[2]] / 2

        J_top = _surface_gas_flux(
            material,
            ξ_top,
            T_cell_top,
            mass_bc_top,
            Δz_half_top,
            t,
            has_darcy_flow(material),
        )
        J_top_thermal = use_transpiration_bc ? J_top : ntuple(_ -> 0.0, Val(NC))

        state_top = solve_boundary_temperature(
            T_cell_top,
            k_top,
            Δz_half_top,
            thermal_bc_top,
            material,
            J_top_thermal,
            T_ref_top,
            t,
            ε_top;
            radiation_model = radiation_model,
            use_blowing = use_blowing_correction,
            J_blowing = J_top,
        )
        solution.surface_T[i] = state_top.T_surface

        # Back (bottom, boundary_id=2) — vents through the user-supplied bottom
        # mass BC (zero for the default impermeable holder; the RHS applies
        # a permeable back-face BC, so the MLR must count it).
        T_cell_bottom = T_profile[1]
        ξ_bottom = ntuple(k -> ξ_profile[k, 1], NC)
        k_bottom = effective_conductivity(ξ_bottom, T_cell_bottom, material)
        ε_bottom = effective_emissivity(ξ_bottom, T_cell_bottom, material)
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
        J_bottom_thermal = use_transpiration_bc ? J_bottom : ntuple(_ -> 0.0, Val(NC))
        state_bottom = solve_boundary_temperature(
            T_cell_bottom,
            k_bottom,
            Δz_half_bottom,
            thermal_bc_bottom,
            material,
            J_bottom_thermal,
            T_ref_bottom,
            t,
            ε_bottom;
            radiation_model = radiation_model,
            use_blowing = use_blowing_correction,
            J_blowing = J_bottom,
        )
        solution.extras.back_temperature[i] = state_bottom.T_surface

        # Total mass (integrate ξ over domain)
        total_m = 0.0
        for j = 1:n_cells_val
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
        ξ_surface = ξ_top

        # Store effective properties
        solution.extras.surface_emissivity[i] = ε_top
        solution.extras.surface_conductivity[i] = k_top

        # Compute heat flux components using SOLVED surface temperature.
        # `blowing_mcp` keeps the convective decomposition consistent with the
        # blowing-corrected flux the Newton applied (zero when blowing off).
        heat_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_top,
            T_surface_solved,
            t,
            ε_top,
            radiation_model,
            state_top.blowing_mcp,
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

        # Transpiration: sensible heat carried out by escaping gas. Closes
        # the surface energy ledger as q_net − q_transpire = q_conduction. Read
        # from the solved boundary state so the diagnostic matches the actual
        # boundary solve (incl. inflow upwinding); zero when transpiration off.
        solution.extras.surface_heat_flux_transpiration[i] = state_top.q_transpired

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
            for j = 1:n_cells_val
                ξ_j = ntuple(k -> ξ_profile[k, j], NC)
                τ_tot += effective_absorption(ξ_j, material) * Geometry.cell_size(mesh, j)
            end
            solution.extras.volumetric_absorption[i] = I_eff * (1.0 - exp(-τ_tot))
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
        ξ_back = ξ_bottom
        back_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_bottom,
            T_back_solved,
            t,
            ε_bottom,
            radiation_model,
            state_bottom.blowing_mcp,
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
        for j = 1:n_cells_val
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

    # For non-ALE mode, thickness is constant (initial mesh thickness)
    initial_thickness = Geometry.domain_length(mesh)
    fill!(solution.extras.thickness_history, initial_thickness)

    return solution
end

"""
	extract_solution_ale!(solution::PyrolysisSolution, raw_solution, problem, layout::StateLayout)

Extract results from ALE ODE solution into PyrolysisSolution structure.

For ALE solutions, node positions are part of the state vector, so we can
reconstruct the geometry at each saved time and compute volumes accordingly.
This ensures proper mass accounting with the evolving mesh.

Note: When surface cell depletion handling is enabled, the state vector length
may vary between saved times. This function detects the current number of cells
from the state vector length and adjusts accordingly.

# Keyword Arguments
- `requested_times::Union{Nothing,Vector{Float64}}`: If provided, filter the raw solution
  to only include times closest to these requested times. This is used when
  depletion handling saves every step and we want specific output times.
	- `radiation_model::RadiationModel`: Radiation model for boundary flux diagnostics.
	- `use_transpiration_bc::Bool`: If false, report boundary temperatures/heat fluxes
	  using the same no-transpiration surface energy balance as the live solve.
	- `use_blowing_correction::Bool`: Apply the same Couette-film blowing rescaling
	  to convective BCs in the diagnostic boundary solve as the live solve did.
	"""
function extract_solution_ale!(
    solution::PyrolysisSolution,
    raw_solution,
    problem,
    layout_initial::StateLayout{NC};
    requested_times::Union{Nothing,AbstractVector{<:Real}} = nothing,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    use_transpiration_bc::Bool = false,
    use_blowing_correction::Bool = false,
) where {NC}
    mesh = problem.mesh
    material = problem.material
    NR = n_reactions(material)

    # Store raw solution for interpolation
    solution.raw_solution = raw_solution

    # Get return code
    solution.retcode = Symbol(raw_solution.retcode)

    # Determine which time indices to extract
    all_times = raw_solution.t
    if requested_times !== nothing && !isempty(requested_times)
        # Find closest actual time for each requested time
        indices_to_extract = Int[]
        for t_req in requested_times
            if t_req > all_times[end]
                break  # Skip requested times beyond simulation end
            end
            idx = argmin(abs.(all_times .- t_req))
            if isempty(indices_to_extract) || idx != indices_to_extract[end]
                push!(indices_to_extract, idx)
            end
        end
        # Always include the final time
        if isempty(indices_to_extract) || indices_to_extract[end] != length(all_times)
            push!(indices_to_extract, length(all_times))
        end
    else
        indices_to_extract = collect(1:length(all_times))
    end

    # Extract times at selected indices
    times = all_times[indices_to_extract]
    solution.t = copy(times)
    n_times = length(times)

    # Pre-allocate output arrays. Under depletion the per-step n_cells/n_nodes
    # shrink along the trajectory; allocate at the MAXIMUM cell count across the
    # extracted snapshots and pad shrunken steps with NaN in T/ξ/z below. Note we
    # cannot rely on `layout_initial.n_cells` here: cell-merge callbacks mutate the
    # workspace layout in place, so by the time extraction runs it may already
    # report the final (merged-down) size. Infer the true maximum from the raw
    # state-vector lengths instead (state layout: [T, ξ, z] plus optional χ).
    has_chi_alloc = !isempty(chi_range(layout_initial))
    per_cell_alloc = NC + 2 + (has_chi_alloc ? 1 : 0)
    n_cells_max = layout_initial.n_cells
    for raw_idx in indices_to_extract
        nc = (length(raw_solution.u[raw_idx]) - 1) ÷ per_cell_alloc
        n_cells_max = max(n_cells_max, nc)
    end
    n_nodes_max = n_cells_max + 1
    solution.T = fill(NaN, n_cells_max, n_times)
    solution.ξ = fill(NaN, NC, n_cells_max, n_times)
    solution.z = fill(NaN, n_nodes_max, n_times)
    solution.surface_T = zeros(n_times)
    solution.mass_loss_rate = zeros(n_times)
    solution.energy_residual = zeros(n_times)
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

    # Extract data at each saved time
    for (out_idx, raw_idx) in enumerate(indices_to_extract)
        t = all_times[raw_idx]
        u = raw_solution.u[raw_idx]

        # Determine current number of cells from state vector length.
        # Layout: [T, ξ, z] plus optional χ.
        state_len = length(u)
        has_chi = !isempty(chi_range(layout_initial))
        n_cells = (state_len - 1) ÷ (NC + 2 + (has_chi ? 1 : 0))
        n_nodes = n_cells + 1

        # Create layout for current size
        layout = StateLayout(NC, n_cells, n_nodes, layout_initial.mode)

        # Extract temperatures (write into T[1:n_cells, out_idx]; remaining
        # rows already NaN from pre-allocation when depletion has shrunk the mesh).
        T_profile = [get_temperature(u, layout, j) for j = 1:n_cells]
        @inbounds for j = 1:n_cells
            solution.T[j, out_idx] = T_profile[j]
        end

        # Extract concentrations
        ξ_profile = zeros(NC, n_cells)
        @inbounds for j = 1:n_cells
            for k = 1:NC
                ξ_profile[k, j] = get_concentration(u, layout, j, k)
                solution.ξ[k, j, out_idx] = ξ_profile[k, j]
            end
        end

        # Extract node positions and compute geometry
        z = get_node_positions(u, layout)
        @inbounds for j = 1:length(z)
            solution.z[j, out_idx] = z[j]
        end
        thickness = z[end] - z[1]
        solution.extras.thickness_history[out_idx] = thickness

        # Compute cell volumes from current node positions and the state-derived A(t).
        A_t = cross_section_area_from_state(mesh, material, u, layout)
        volumes = compute_volumes_from_state(u, layout, A_t)

        # ================================================================
        # Solve for true boundary temperatures. We reconstruct the per-species
        # gas mass flux at the top face the same way the live mass BC does,
        # so the boundary-temperature solve here matches the integrator's
        # transpiration BC. The back face has no mass BC (insulated holder),
        # so J_bottom is zero by construction.
        # ================================================================
        T_ref_top = ambient_temperature(thermal_bc_top, t, Float64(T_REF))
        T_ref_bottom = ambient_temperature(thermal_bc_bottom, t, Float64(T_REF))

        # Surface (top, boundary_id=1) — geometry, properties, gas flux
        T_cell_top = T_profile[n_cells]
        ξ_top = ntuple(k -> ξ_profile[k, n_cells], NC)
        k_top = effective_conductivity(ξ_top, T_cell_top, material)
        ε_top = effective_emissivity(ξ_top, T_cell_top, material)
        Δz_half_top = (z[n_cells+1] - z[n_cells]) / 2

        J_top = _surface_gas_flux(
            material,
            ξ_top,
            T_cell_top,
            mass_bc_top,
            Δz_half_top,
            t,
            has_darcy_flow(material),
        )
        J_top_thermal = use_transpiration_bc ? J_top : ntuple(_ -> 0.0, Val(NC))

        state_top = solve_boundary_temperature(
            T_cell_top,
            k_top,
            Δz_half_top,
            thermal_bc_top,
            material,
            J_top_thermal,
            T_ref_top,
            t,
            ε_top;
            radiation_model = radiation_model,
            use_blowing = use_blowing_correction,
            J_blowing = J_top,
        )
        solution.surface_T[out_idx] = state_top.T_surface

        # Back (bottom, boundary_id=2) — vents through the user-supplied bottom
        # mass BC (zero for the default impermeable holder; the RHS applies
        # a permeable back-face BC, so the MLR must count it).
        T_cell_bottom = T_profile[1]
        ξ_bottom = ntuple(k -> ξ_profile[k, 1], NC)
        k_bottom = effective_conductivity(ξ_bottom, T_cell_bottom, material)
        ε_bottom = effective_emissivity(ξ_bottom, T_cell_bottom, material)
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
        J_bottom_thermal = use_transpiration_bc ? J_bottom : ntuple(_ -> 0.0, Val(NC))
        state_bottom = solve_boundary_temperature(
            T_cell_bottom,
            k_bottom,
            Δz_half_bottom,
            thermal_bc_bottom,
            material,
            J_bottom_thermal,
            T_ref_bottom,
            t,
            ε_bottom;
            radiation_model = radiation_model,
            use_blowing = use_blowing_correction,
            J_blowing = J_bottom,
        )
        solution.extras.back_temperature[out_idx] = state_bottom.T_surface

        # Total mass (integrate ξ × V using current geometry)
        total_m = 0.0
        for j = 1:n_cells
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
        ξ_surface = ξ_top

        # Store effective properties
        solution.extras.surface_emissivity[out_idx] = ε_top
        solution.extras.surface_conductivity[out_idx] = k_top

        # Compute heat flux components using SOLVED surface temperature.
        # `blowing_mcp` keeps the convective decomposition consistent with the
        # blowing-corrected flux the Newton applied (zero when blowing off).
        heat_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_top,
            T_surface_solved,
            t,
            ε_top,
            radiation_model,
            state_top.blowing_mcp,
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

        # Transpiration: sensible heat carried out by escaping gas. Closes
        # the surface energy ledger as q_net − q_transpire = q_conduction. Read
        # from the solved boundary state so the diagnostic matches the actual
        # boundary solve (incl. inflow upwinding); zero when transpiration off.
        solution.extras.surface_heat_flux_transpiration[out_idx] = state_top.q_transpired

        # Volumetric radiative absorption (BEER_LAMBERT only) on the live ALE
        # geometry: Σ q'''·Δz = α_surf·q_ext·(1 − e^{−τ_tot}); see the
        # Eulerian extractor for the convention.
        if radiation_model == BEER_LAMBERT
            I_inc, has_explicit_α = extract_radiative_flux_from_bcs(bc_set, Val(:top), t)
            I_eff = has_explicit_α ? I_inc : ε_top * I_inc
            τ_tot = 0.0
            for j = 1:n_cells
                ξ_j = ntuple(k -> ξ_profile[k, j], NC)
                τ_tot += effective_absorption(ξ_j, material) * (z[j+1] - z[j])
            end
            solution.extras.volumetric_absorption[out_idx] = I_eff * (1.0 - exp(-τ_tot))
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
        ξ_back = ξ_bottom
        back_fluxes = compute_boundary_heat_fluxes(
            thermal_bc_bottom,
            T_back_solved,
            t,
            ε_bottom,
            radiation_model,
            state_bottom.blowing_mcp,
        )
        solution.extras.back_heat_flux[out_idx] =
            thermal_bc_bottom isa DirichletThermalBC ? state_bottom.q_conducted :
            back_fluxes.net

        # Gas generation rate: compute from reaction rates (physics-based)
        # This is NOT the actual mass loss rate (surface flux), but rather the
        # volumetric production of gases throughout the sample.
        gen_rate = 0.0
        for j = 1:n_cells
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
