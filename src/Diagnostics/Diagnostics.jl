"""
    Diagnostics

Conservation diagnostics.

- [`MassConservationTracker`](@ref) — per-component mass ledger.
- [`EnergyConservationTracker`](@ref) — sensible-energy ledger of the
  temperature equation (`E = Σ ρcp·T·V` over every component), balanced
  against the terms the discretized RHS assembles.
- [`DiagnosticsState`](@ref) and [`record_substep!`](@ref) — advance both
  ledgers by one step's boundary fluxes and volumetric sources and log a row
  every `stride` steps. The test suite drives them step by step to check
  conservation; `solve` does not.
- [`DiagnosticsCallback`](@ref) — `DiscreteCallback` builder that runs a
  `record!(u, t)` closure on every accepted step. `solve` uses it to log
  [`snapshot_diagnostics_row`](@ref)s (the sensible energy at the step) when
  diagnostics are requested.
"""
module Diagnostics

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties

include("conservation_checks.jl")
include("callback.jl")

# Exports — trackers
export MassConservationTracker, EnergyConservationTracker
export compute_domain_mass, compute_sensible_energy
export update_mass_tracker!, update_energy_tracker!
export check_mass_conservation, check_energy_conservation
export mass_conservation_report, energy_conservation_report

# Exports — callback infrastructure
export DiagnosticsState,
    DiagnosticsRow, DiagnosticsCallback, record_substep!, snapshot_diagnostics_row

end # module Diagnostics
