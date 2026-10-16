"""
Diagnostics callback infrastructure: the typed log row, the ledger state
advanced by [`record_substep!`](@ref), the snapshot row `solve` logs, and the
`DiscreteCallback` builder.
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

Ledger state advanced by [`record_substep!`](@ref). Holds:
- the running mass and energy trackers,
- a row log written every `stride` recorded steps,
- the time of the last recorded step (for the substep `dt`).
"""
mutable struct DiagnosticsState{NC}
    mass::MassConservationTracker{NC}
    energy::EnergyConservationTracker
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
    return DiagnosticsState{NC}(
        mass_tracker,
        energy_tracker,
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
    state::DiagnosticsState{NC},
    mesh::Geometry.Mesh1D{NC},
    material;
    q_cond_boundary::Float64,
    Q_rxn_total::Float64,
    S_conv_total::Float64,
    J_boundary::Matrix{Float64},
    A::Float64,
    t::Float64,
) where {NC}

    dt = t - state.last_t
    state.last_t = t

    if dt <= 0
        return nothing  # callback fired at duplicate t (e.g. the initial point)
    end

    state.step_counter += 1

    # Update the mass and energy trackers.
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

Snapshot at a single point: the sensible energy. The ledger fields
(cumulative residual, mass error) are `NaN`; only [`record_substep!`](@ref)
fills them. This is the row `solve` logs.
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
log when appropriate; `solve` pushes [`snapshot_diagnostics_row`](@ref)s.

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
