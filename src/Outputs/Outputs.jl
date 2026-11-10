"""
	Outputs

The solution: the saved times and, at each, the
fields reconstructed from the state through the state mapping and the
kernels; its accessors; and the conservation ledger of the live region with
its per-step recorder and the audit that closes it.
"""
module Outputs

using ..Meshes
using ..States
using ..States: state_length
using ..Materials: is_gas, T_REF, gas_component_indices
using ..Physics: pore_state
using ..BoundaryConditions: compute_boundary_heat_fluxes, BlowingFace, DirichletThermalBC
using ..Kernels
using ..Assembly
using ..Deactivation

include("gas_diagnostics.jl")
include("solution.jl")
include("conservation.jl")
include("accessors.jl")

export GasDiagnostics, get_gas_composition, get_boundary_exchange, gas_species_audit
export PyrolysisSolution, SolutionExtras, SolverStatistics
export succeeded, completed
export get_profile, interior_events, displaced_residue
export ConservationRecorder, conservation_audit, ledger_rates, live_inventories

end # module Outputs
