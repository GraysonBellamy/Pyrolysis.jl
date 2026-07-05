"""
    Diagnostics

First-class runtime conservation diagnostics for the Pyrolysis.jl
energy formulation.

Two trackers + a callback:

- [`MassConservationTracker`](@ref) — per-component mass ledger.
- [`EnergyConservationTracker`](@ref) — both the **discrete-form** residual
  (closes to integrator tolerance by construction) and the **physical-form**
  residual (quantifies the cost of the quasi-steady-gas approximation that
  Formulation A makes by excluding gases from `ρcp_eff`).
- [`ReactionMassFlowTracker`](@ref) — per-reaction mass balance.
- [`DiagnosticsCallback`](@ref) — `DiscreteCallback` builder that fires every
  `stride` steps and at the final time, accumulating the substep contributions
  into the trackers and emitting a row to a `Vector{NamedTuple}` log.

The callback is the visibility mechanism that makes the FDS-convention
approximations (gas-excluded `ρcp_eff`, dropped `T·cp·∇·J` term) auditable
on every solve. The log is connected to `PyrolysisSolution.extras`.
"""
module Diagnostics

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties

include("conservation_checks.jl")
include("callback.jl")

# Exports — trackers
export MassConservationTracker, EnergyConservationTracker, ReactionMassFlowTracker
export compute_domain_mass, compute_sensible_energy, compute_total_sensible_energy
export update_mass_tracker!, update_energy_tracker!, update_reaction_tracker!
export check_mass_conservation, check_energy_conservation, check_reaction_mass_balance
export mass_conservation_report, energy_conservation_report, reaction_mass_flow_report

# Exports — callback infrastructure
export DiagnosticsState,
    DiagnosticsRow, DiagnosticsCallback, record_substep!, snapshot_diagnostics_row

end # module Diagnostics
