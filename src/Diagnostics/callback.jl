"""
Diagnostics callback infrastructure.

Builds a `DiscreteCallback` that fires every
`stride` accepted timesteps (and at the final time) to:

1. Update the running mass / energy / reaction trackers with this substep's
   contributions, and
2. Emit a typed `NamedTuple` row to `log::Vector{NamedTuple}` capturing the
   energy and mass ledger residuals at this point in the simulation.

The log is connected to `PyrolysisSolution.extras`.
"""

using SciMLBase: DiscreteCallback, derivative_discontinuity!
using FunctionWrappers: FunctionWrapper

# Row type — concretely typed so consumers don't have to deal with
# Dict{Symbol,Any}. Adding a new diagnostic adds a field here, not a key.
const DiagnosticsRow = NamedTuple{
    (:t, :E_sensible, :energy_residual, :energy_relative, :mass_max_relative_error),
    Tuple{Float64,Float64,Float64,Float64,Float64},
}

"""
    DiagnosticsState

Mutable container threaded through the callback. Holds:
- the running trackers,
- a row log written every `stride` accepted steps,
- the last time at which the callback ran (for substep dt computation).

Built once at the start of `solve()` and lives on the `Workspace` (or as a
free-standing object passed via `cb_kwargs`).
"""
mutable struct DiagnosticsState{NC,NR}
    mass::MassConservationTracker{NC}
    energy::EnergyConservationTracker
    reactions::ReactionMassFlowTracker{NC,NR}
    stride::Int
    step_counter::Int
    last_t::Float64
    log::Vector{DiagnosticsRow}
end

"""
    DiagnosticsState(mesh, material; stride=10) -> DiagnosticsState

Initialize the diagnostics state with fresh trackers seeded from the current
mesh state (this becomes `t=0` baseline).
"""
function DiagnosticsState(mesh::Geometry.Mesh1D{NC}, material; stride::Int = 10) where {NC}
    mass_tracker = MassConservationTracker(mesh, material)
    energy_tracker = EnergyConservationTracker(mesh, material)
    reaction_tracker = ReactionMassFlowTracker(material)
    NR = length(material.reactions)
    return DiagnosticsState{NC,NR}(
        mass_tracker,
        energy_tracker,
        reaction_tracker,
        stride,
        0,
        0.0,
        DiagnosticsRow[],
    )
end

"""
    record_substep!(state::DiagnosticsState, mesh, material;
                    q_cond_boundary, Q_rxn_total, S_conv_total,
                    J_boundary, A, t)

Update the trackers with one substep's contributions and append a typed row
to the log if `step_counter % stride == 0`.

Substep size `dt = t − state.last_t` is computed automatically.
"""
function record_substep!(
    state::DiagnosticsState{NC,NR},
    mesh::Geometry.Mesh1D{NC},
    material;
    q_cond_boundary::Float64,
    Q_rxn_total::Float64,
    S_conv_total::Float64,
    J_boundary::Matrix{Float64},
    A::Float64,
    t::Float64,
) where {NC,NR}

    dt = t - state.last_t
    state.last_t = t

    if dt <= 0
        return nothing  # callback fired at duplicate t (e.g. the initial point)
    end

    state.step_counter += 1

    # Update the mass, energy, and reaction trackers.
    update_mass_tracker!(state.mass, J_boundary, A, dt)
    update_energy_tracker!(
        state.energy;
        q_cond_boundary = q_cond_boundary,
        Q_rxn_total = Q_rxn_total,
        S_conv_total = S_conv_total,
        A = A,
        dt = dt,
    )

    # Emit a row only on stride boundaries.
    if state.step_counter % state.stride == 0
        ec = check_energy_conservation(state.energy, mesh, material)
        mc = check_mass_conservation(state.mass, mesh, material)
        row = DiagnosticsRow((
            t = t,
            E_sensible = ec.current_energy,
            energy_residual = ec.residual,
            energy_relative = ec.relative,
            mass_max_relative_error = mc.max_error,
        ))
        push!(state.log, row)
    end
    return nothing
end

"""
    snapshot_diagnostics_row(t, mesh, material) -> DiagnosticsRow

Cumulative-free snapshot at a single point: the sensible energy. The ledger
fields (cumulative residual, mass error) are filled with `NaN` to mark them
as not-yet-computed. This is what the snapshot-mode callback writes; the
cumulative form requires per-substep accumulation that the unified residual
layer wires through `record_substep!`.
"""
function snapshot_diagnostics_row(t::Float64, mesh::Geometry.Mesh1D, material)
    return DiagnosticsRow((
        t = t,
        E_sensible = compute_sensible_energy(mesh, material),
        energy_residual = NaN,
        energy_relative = NaN,
        mass_max_relative_error = NaN,
    ))
end

"""
    DiagnosticsCallback(record!)

Construct a `DiscreteCallback` that fires every accepted step and runs
`record!(u, t)` on the integrator's state `u::Vector{Float64}` and time
`t::Float64`. It only observes the state, so it clears the
derivative-discontinuity flag. Call sites supply a `record!`
closure that has access to the live workspace and pushes a row to the
log when appropriate. This module ships [`snapshot_diagnostics_row`](@ref)
as the simplest builder; a full ledger callback that calls
[`record_substep!`](@ref) requires per-substep accumulators that the unified
residual layer threads through.

`record!` is held behind a fixed-signature `FunctionWrapper`, so the
callback's type does not depend on the closure's and the integrator it joins
keeps a problem-independent type.
"""
function DiagnosticsCallback(record!)
    record_wrapped! = FunctionWrapper{Nothing,Tuple{Vector{Float64},Float64}}(record!)
    condition(_, _, _) = true
    function affect!(integrator)
        record_wrapped!(integrator.u, integrator.t)
        derivative_discontinuity!(integrator, false)
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end
