"""
	Solver

Solution algorithms for the Pyrolysis.jl package.

This module provides:
- Main solve interface
- Operator splitting
- Callbacks
- Steady-state solver

ALE mesh motion is integrated into the unified residual / workspace pair.
This ensures proper mass conservation.
"""
module Solver

using ..Materials
using ..BoundaryConditions
using ..Materials: n_components
using ..Geometry
using ..Properties
using ..Physics
using ..Adaptivity
using ..Discretization
using ..Jacobian
using ..Jacobian:
    JacobianSpec,
    JacobianBackend,
    compute_jacobian!,
    Structured,
    DirectSolve,
    ApproxSolve,
    StructuredSolve,
    DenseAD_ForwardDiff,
    verify_jacobian_at,
    sparsity_profile,
    materializes_geometry_block
using ..Problem
using ..Diagnostics: DiagnosticsRow, DiagnosticsCallback, snapshot_diagnostics_row
using OrdinaryDiffEq
using OrdinaryDiffEq: CallbackSet
# The default integrator KenCarp4 lives in the SDIRK sublibrary, outside the
# umbrella package's default solver set.
using OrdinaryDiffEqSDIRK: KenCarp4
using SciMLBase: terminate!, derivative_discontinuity!, set_proposed_dt!, get_proposed_dt
using SciMLBase: FullSpecialize
using SciMLBase: DiscreteCallback, ContinuousCallback
using Sparspak  # Loads LinearSolve's Sparspak extension, so `SparspakFactorization()` is selectable
using ADTypes: AutoFiniteDiff
using DiffEqCallbacks: PeriodicCallback
using ProgressMeter: Progress, update!, finish!

# Include subfiles
include("callbacks.jl")
include("solve.jl")
include("ilu_gmres.jl")
include("parameter_sweep.jl")

# Exports
# Main solve interface
export solve
export ILUGMRESFactorization, ILUBiCGSTABFactorization
export parameter_sweep

# Solution utilities
export get_profile, get_cell_positions, interpolate_solution

# Callbacks
export create_thickness_termination_callback
export create_depletion_callback, DepletionCallbackState, DT_AFTER_PARK, n_live_cells
export create_step_update_callback
export create_newton_failure_cap_callback, FAILURE_MEMORY_GROWTH, FAILURE_MEMORY_STEPS

# Thickness-based detection (used by depletion callback)
export is_cell_thin, get_cell_thickness, find_thin_cells

end # module Solver
