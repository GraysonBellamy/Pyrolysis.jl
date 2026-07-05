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
    verify_jacobian_at
using ..Problem
using ..Diagnostics: DiagnosticsRow, DiagnosticsCallback, snapshot_diagnostics_row
using OrdinaryDiffEq
using OrdinaryDiffEq: CallbackSet
# OrdinaryDiffEq v7 split solver families into sublibraries: only the default set
# (Tsit5/Rodas5P/FBDF/…) stays in the umbrella, so the default integrator KenCarp4
# (SDIRK) must be imported explicitly. See OrdinaryDiffEq NEWS.md (v7).
using OrdinaryDiffEqSDIRK: KenCarp4
using SciMLBase: terminate!, u_modified!  # owned by SciMLBase; OrdinaryDiffEq v7 no longer re-exports them
using Sparspak  # Required for SparspakFactorization (default linear solver)
using ADTypes: AutoFiniteDiff
using DiffEqCallbacks: PeriodicCallback, ContinuousCallback, DiscreteCallback
using ProgressMeter: Progress, update!, finish!

# Include subfiles
include("solve.jl")
include("ilu_gmres.jl")
include("callbacks.jl")
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
export create_depletion_callback, DepletionCallbackState
export get_depletion_callback_stats
export create_step_update_callback

# Thickness-based detection (used by depletion callback)
export is_cell_thin, get_cell_thickness, find_thin_cells

end # module Solver
