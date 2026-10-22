"""
	Solver

Solution algorithms for the Pyrolysis.jl package.

This module provides:
- the `solve` interface and solution extraction;
- the integrator callbacks (step update, cell parking, Newton-failure cap,
  diagnostics);
- the ILU-GMRES linear solver and `parameter_sweep`.

ALE mesh motion is part of the residual's state (node positions), so mesh and
material move together in one implicit system.
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
    compute_jacobian!,
    Structured,
    DirectSolve,
    ApproxSolve,
    DenseAD_ForwardDiff,
    verify_jacobian_at,
    sparsity_profile,
    materializes_geometry_block
using ..Problem
# The accessors and the sweep are the finite-volume core's; the old solution
# and problem add their methods.
import ..Outputs: get_profile, get_cell_positions, interpolate_solution
import ..Integration: parameter_sweep
using ..Diagnostics: DiagnosticsRow, DiagnosticsCallback, snapshot_diagnostics_row
using OrdinaryDiffEq
using OrdinaryDiffEq: CallbackSet
# The default integrator KenCarp4 lives in the SDIRK sublibrary, outside the
# umbrella package's default solver set.
using OrdinaryDiffEqSDIRK: KenCarp4
using SciMLBase: terminate!, derivative_discontinuity!, set_proposed_dt!, get_proposed_dt
using SciMLBase: FullSpecialize, NullParameters
using SciMLBase: DiscreteCallback, ContinuousCallback
using Sparspak  # Loads LinearSolve's Sparspak extension, so `SparspakFactorization()` is selectable
using ADTypes: AutoFiniteDiff
using ProgressMeter: Progress, update!, finish!
using FunctionWrappers: FunctionWrapper

# Include subfiles
include("integrator_functions.jl")
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
export is_cell_thin, find_thin_cells

end # module Solver
