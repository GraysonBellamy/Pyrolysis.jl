"""
	Solvers

Solution algorithms for the Pyrolysis.jl package.

This module provides:
- Main solve interface
- Operator splitting
- Callbacks
- Steady-state solver

ALE mesh motion is now integrated into the ODE function (PyrolysisODEFunctionALE)
rather than using callbacks. This ensures proper mass conservation.
"""
module Solvers

using ..Core
using ..Core: n_components
using ..Geometry
using ..Properties
using ..Physics
using ..ALE
using ..Adaptivity
using ..Discretization
using ..TimeIntegration
using ..Problem
using OrdinaryDiffEq
using OrdinaryDiffEq: CallbackSet, terminate!
using Sparspak  # Required for SparspakFactorization (default linear solver)
using ADTypes: AutoFiniteDiff
using DiffEqCallbacks: PeriodicCallback, ContinuousCallback
using ProgressMeter: Progress, update!, finish!

# Include subfiles
include("solve.jl")
include("ilu_gmres.jl")
include("callbacks.jl")

# Exports
# Main solve interface
export solve
export ILUGMRESFactorization, ILUBiCGSTABFactorization

# Solution utilities
export get_profile, get_cell_positions, interpolate_solution

# Callbacks
export create_thickness_termination_callback
export create_depletion_callback, DepletionCallbackState
export get_depletion_callback_stats

# Thickness-based detection (used by depletion callback)
export is_cell_thin, get_cell_thickness, find_thin_cells

end # module Solvers
