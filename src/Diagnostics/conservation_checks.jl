"""
Conservation checks for mass and energy.

This module provides runtime verification of conservation laws to detect
numerical errors and validate physics implementations.
"""

using StaticArrays
using Printf

# ============================================================================
# Mass Conservation Tracking
# ============================================================================

"""
    MassConservationTracker{NC}

Tracks cumulative mass for each component throughout simulation.

# Fields
- `initial_mass::NTuple{NC, Float64}`: Initial mass per component [kg]
- `cumulative_efflux::NTuple{NC, Float64}`: Cumulative mass leaving domain [kg]
- `last_time::Float64`: Time of last update [s]
"""
mutable struct MassConservationTracker{NC}
    initial_mass::NTuple{NC,Float64}
    cumulative_efflux::NTuple{NC,Float64}
    last_time::Float64
end

"""
    MassConservationTracker(mesh::Mesh1D{NC}, material) where NC

Initialize tracker with current domain mass.
"""
function MassConservationTracker(mesh::Geometry.Mesh1D{NC}, material) where {NC}
    initial_mass = compute_domain_mass(mesh, material)
    cumulative_efflux = ntuple(_ -> 0.0, Val(NC))
    return MassConservationTracker{NC}(initial_mass, cumulative_efflux, 0.0)
end

"""
    compute_domain_mass(mesh::Mesh1D{NC}, material) -> NTuple{NC, Float64}

Compute total mass of each component in the domain.

Mass for component j: mⱼ = Σᵢ ξⱼ(i) × V(i)

Only includes active cells (excludes refined-away cells in AMR).
"""
function compute_domain_mass(mesh::Geometry.Mesh1D{NC}, material) where {NC}
    mass = zeros(MVector{NC,Float64})

    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end

        V = mesh.cell_volumes[i]
        ξ = mesh.cell_states[i].ξ

        for j = 1:NC
            mass[j] += ξ[j] * V
        end
    end

    return NTuple{NC,Float64}(mass)
end

"""
    update_mass_tracker!(tracker::MassConservationTracker{NC},
                         J_boundary::Matrix{Float64},
                         A::Float64, dt::Float64) where NC

Update tracker with boundary mass fluxes.

# Arguments
- `tracker`: The mass conservation tracker
- `J_boundary`: Boundary mass fluxes [NC × 2] for top/bottom boundaries [kg/(m²·s)]
- `A`: Cross-sectional area [m²]. Under variable cross-section,
  pass the trapezoidal-in-A average over the substep,
  `A = 0.5 * (A_t_pre + A_t_post)`. With identity law / non-ALE, `A` is
  constant and equals `mesh.cross_section_area_initial`.
- `dt`: Time step [s]

# Notes
- Convention: Positive J = flux into domain (inflow)
- Mass leaving = -J × A × dt (negative J = outflow)
- Typical case: J_top < 0 (gas escaping), J_bottom = 0 (impermeable)
"""
function update_mass_tracker!(
    tracker::MassConservationTracker{NC},
    J_boundary::Matrix{Float64},
    A::Float64,
    dt::Float64,
) where {NC}
    new_efflux = MVector{NC,Float64}(tracker.cumulative_efflux)

    @inbounds for j = 1:NC
        # Sum fluxes from top (index 1) and bottom (index 2) boundaries
        J_top = J_boundary[j, 1]
        J_bottom = J_boundary[j, 2]

        # Mass leaving domain = -J × A × dt (positive J = into domain)
        mass_out = -(J_top + J_bottom) * A * dt
        new_efflux[j] += mass_out
    end

    tracker.cumulative_efflux = NTuple{NC,Float64}(new_efflux)
    tracker.last_time += dt

    return nothing
end

"""
    check_mass_conservation(tracker::MassConservationTracker{NC},
                            mesh::Mesh1D{NC}, material;
                            tol::Float64=1e-10) -> NamedTuple

Check mass conservation and return detailed diagnostics.

Conservation equation:
    m_current(t) + m_escaped(t) = m_initial

Or equivalently:
    m_current(t) = m_initial - m_escaped(t)

# Arguments
- `tracker`: Mass conservation tracker
- `mesh`: Current mesh state
- `material`: Material definition
- `tol`: Relative error tolerance (default: 1e-10)

# ALE note
Every component rides the mesh (its only mesh-motion term is the dilation
`−ξθ`, which conserves the cell mass `ξV`), so under ALE the mass of every
component changes only through the transport flux `J` at the boundary faces
— exactly what `expected_mass` integrates. A fully closed (impermeable)
column conserves each component to integrator tolerance (V-CPL-03 P4).

# Returns
NamedTuple with:
- `conserved::Bool`: True if mass is conserved within tolerance
- `current_mass::NTuple{NC, Float64}`: Current mass in domain [kg]
- `expected_mass::NTuple{NC, Float64}`: Expected mass (initial - efflux) [kg]
- `relative_error::NTuple{NC, Float64}`: Relative error per component
- `max_error::Float64`: Maximum relative error across components
"""
function check_mass_conservation(
    tracker::MassConservationTracker{NC},
    mesh::Geometry.Mesh1D{NC},
    material;
    tol::Float64 = 1e-10,
) where {NC}
    current_mass = compute_domain_mass(mesh, material)

    # Expected mass = initial - cumulative efflux
    expected_mass =
        ntuple(j -> tracker.initial_mass[j] - tracker.cumulative_efflux[j], Val(NC))

    # Compute relative errors
    rel_error = ntuple(Val(NC)) do j
        if abs(tracker.initial_mass[j]) > eps()
            abs(current_mass[j] - expected_mass[j]) / tracker.initial_mass[j]
        else
            # For components with zero initial mass, use absolute error
            abs(current_mass[j] - expected_mass[j])
        end
    end

    max_error = maximum(rel_error)
    conserved = max_error < tol

    return (
        conserved = conserved,
        current_mass = current_mass,
        expected_mass = expected_mass,
        relative_error = rel_error,
        max_error = max_error,
    )
end

"""
    mass_conservation_report(tracker::MassConservationTracker{NC},
                             mesh::Mesh1D{NC}, material) -> String

Generate a human-readable mass conservation report.

Useful for post-processing and debugging.
"""
function mass_conservation_report(
    tracker::MassConservationTracker{NC},
    mesh::Geometry.Mesh1D{NC},
    material,
) where {NC}
    result = check_mass_conservation(tracker, mesh, material)

    lines = String[]
    push!(lines, "=" ^ 60)
    push!(lines, "MASS CONSERVATION REPORT")
    push!(lines, "=" ^ 60)
    push!(lines, @sprintf("Time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 60)

    for j = 1:NC
        comp_name = string(material.components[j].name)
        push!(lines, @sprintf("Component %d (%s):", j, comp_name))
        push!(lines, @sprintf("  Initial mass:    %12.6e kg", tracker.initial_mass[j]))
        push!(lines, @sprintf("  Cumulative out:  %12.6e kg", tracker.cumulative_efflux[j]))
        push!(lines, @sprintf("  Expected mass:   %12.6e kg", result.expected_mass[j]))
        push!(lines, @sprintf("  Current mass:    %12.6e kg", result.current_mass[j]))
        push!(lines, @sprintf("  Relative error:  %12.6e", result.relative_error[j]))
        push!(lines, "")
    end

    push!(lines, "-" ^ 60)
    status = result.conserved ? "PASSED" : "FAILED"
    push!(lines, @sprintf("Overall: %s (max error = %.2e)", status, result.max_error))
    push!(lines, "=" ^ 60)

    return join(lines, "\n")
end

# ============================================================================
# Energy Conservation Tracking
# ============================================================================

"""
    EnergyConservationTracker

Sensible-energy ledger of the temperature equation. The stored quantity is

    E(t) = Σ_i ρcp_i · T_i · V_i,    ρcp = Σ_j ξ_j cp_j(T)  (every component)

and the ledger accumulates the terms the discretized RHS actually assembles:

    ΔE(t) = ∫ q_cond_boundary · A · dt
          + ∫ Q_rxn dV · dt
          + ∫ S_conv dV · dt

`q_cond_boundary` is the conducted boundary flux — the energy the boundary
cell receives. The surface energy balance has no gas-enthalpy term, so it
equals the applied boundary flux. The residual closes to roughly
integrator tolerance × |E(0)| for constant composition and volume; with
active reactions or a moving mesh it also carries the composition-change
term `∫ T ∂(ρcp)/∂t dV dt` and the volume term `∫ ρcp T ∂V/∂t dt`, which
the per-substep accumulation below does not reconstruct. NOTE: because it
is built from the terms the RHS *actually* assembles, this ledger closes
even if the RHS itself double-counts a mechanism — it audits discrete
bookkeeping, not the physical model.

# Fields
- `initial_energy::Float64`                 — `Σ ρcp · T · V` at t=0.
- `cumulative_boundary_conduction::Float64` — `∫ q_cond_boundary · A · dt`.
- `cumulative_reaction_heat::Float64`       — `∫ Q_rxn dV · dt`.
- `cumulative_advection_sink::Float64`      — `∫ ∫ S_conv dV · dt`.
- `last_time::Float64`                      — Time of last update [s].
"""
mutable struct EnergyConservationTracker
    initial_energy::Float64
    cumulative_boundary_conduction::Float64
    cumulative_reaction_heat::Float64
    cumulative_advection_sink::Float64
    last_time::Float64
end

"""
    EnergyConservationTracker(mesh::Mesh1D, material)

Initialize tracker with the current sensible energy.
"""
function EnergyConservationTracker(mesh::Geometry.Mesh1D, material)
    E0 = compute_sensible_energy(mesh, material)
    return EnergyConservationTracker(E0, 0.0, 0.0, 0.0, 0.0)
end

"""
    compute_sensible_energy(mesh::Mesh1D{NC}, material) -> Float64

Sensible heat `E = Σ_i ρcp_i · T_i · V_i` with the mixture `ρcp` (every
component, gas included — `Properties.effective_heat_capacity`). Only active
cells (AMR-aware).
"""
function compute_sensible_energy(mesh::Geometry.Mesh1D{NC}, material) where {NC}
    E = 0.0
    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end

        state = mesh.cell_states[i]
        T = state.T
        ξ = state.ξ
        V = mesh.cell_volumes[i]

        ρcp = Properties.effective_heat_capacity(ξ, T, material)
        E += ρcp * T * V
    end
    return E
end

"""
    update_energy_tracker!(tracker; q_cond_boundary, Q_rxn_total,
                           S_conv_total, A, dt)

Update tracker with one substep's contributions.

All `*_total` arguments are **volume-integrals over the cell-source field**
(`Σ_i field_i · V_i`) computed from the active RHS state. `q_cond_boundary`
is a face-flux value [W/m²] (multiplied by A and dt internally).

# Keyword arguments
- `q_cond_boundary::Float64` — net conducted boundary flux [W/m²], positive into slab.
- `Q_rxn_total::Float64`     — `∫ Q_rxn dV` [W].
- `S_conv_total::Float64`    — `∫ S_conv dV` [W] (typically negative).
- `A::Float64`               — Cross-sectional area at this substep [m²].
- `dt::Float64`              — Substep size [s].
"""
function update_energy_tracker!(
    tracker::EnergyConservationTracker;
    q_cond_boundary::Float64,
    Q_rxn_total::Float64,
    S_conv_total::Float64,
    A::Float64,
    dt::Float64,
)
    tracker.cumulative_boundary_conduction += q_cond_boundary * A * dt
    tracker.cumulative_reaction_heat += Q_rxn_total * dt
    tracker.cumulative_advection_sink += S_conv_total * dt
    tracker.last_time += dt
    return nothing
end

"""
    check_energy_conservation(tracker, mesh, material; tol=1e-6) -> NamedTuple

Compute the energy-ledger residual `E(t) − E(0) − Σ contributions`.

# Returns
NamedTuple with:
- `residual::Float64`        — `E(t) − E(0) − Σ contributions` [J]. `≈ 0` to integrator tolerance for constant composition and volume.
- `relative::Float64`        — relative magnitude vs `|E(0)|`.
- `conserved::Bool`          — `relative < tol`.
- `current_energy::Float64`  — current sensible heat [J].
"""
function check_energy_conservation(
    tracker::EnergyConservationTracker,
    mesh::Geometry.Mesh1D,
    material;
    tol::Float64 = 1e-6,
)
    E_now = compute_sensible_energy(mesh, material)

    expected =
        tracker.initial_energy +
        tracker.cumulative_boundary_conduction +
        tracker.cumulative_reaction_heat +
        tracker.cumulative_advection_sink
    residual = E_now - expected

    relative = if abs(tracker.initial_energy) > eps()
        abs(residual) / abs(tracker.initial_energy)
    else
        abs(residual)
    end

    return (
        residual = residual,
        relative = relative,
        conserved = relative < tol,
        current_energy = E_now,
    )
end

"""
    energy_conservation_report(tracker::EnergyConservationTracker,
                               mesh::Mesh1D, material) -> String

Generate a human-readable energy conservation report.

Useful for post-processing and debugging.
"""
function energy_conservation_report(
    tracker::EnergyConservationTracker,
    mesh::Geometry.Mesh1D,
    material,
)
    result = check_energy_conservation(tracker, mesh, material)

    lines = String[]
    push!(lines, "=" ^ 70)
    push!(lines, "ENERGY CONSERVATION REPORT")
    push!(lines, "=" ^ 70)
    push!(lines, @sprintf("Time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 70)
    push!(lines, @sprintf("  E(0):                       %12.6e J", tracker.initial_energy))
    push!(lines, @sprintf("  E(t):                       %12.6e J", result.current_energy))
    push!(
        lines,
        @sprintf(
            "  ∫ q_cond_boundary · A · dt: %12.6e J",
            tracker.cumulative_boundary_conduction
        )
    )
    push!(
        lines,
        @sprintf(
            "  ∫ Q_rxn dV · dt:            %12.6e J",
            tracker.cumulative_reaction_heat
        )
    )
    push!(
        lines,
        @sprintf(
            "  ∫∫ S_conv dV · dt:          %12.6e J",
            tracker.cumulative_advection_sink
        )
    )
    push!(
        lines,
        @sprintf(
            "  Residual:                   %12.6e J  (rel. %.2e) %s",
            result.residual,
            result.relative,
            result.conserved ? "✓" : "✗"
        )
    )
    push!(lines, "=" ^ 70)
    return join(lines, "\n")
end
