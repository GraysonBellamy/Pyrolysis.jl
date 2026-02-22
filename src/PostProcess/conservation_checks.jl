"""
Conservation checks for mass and energy.

Reference: Design Document Section 10.3
Physics Audit Recommendations (2026)

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
    initial_mass::NTuple{NC, Float64}
    cumulative_efflux::NTuple{NC, Float64}
    last_time::Float64
end

"""
    MassConservationTracker(mesh::Mesh1D{NC}, material) where NC

Initialize tracker with current domain mass.
"""
function MassConservationTracker(mesh::Geometry.Mesh1D{NC}, material) where NC
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
function compute_domain_mass(mesh::Geometry.Mesh1D{NC}, material) where NC
    mass = zeros(MVector{NC, Float64})

    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end

        V = mesh.cell_volumes[i]
        ξ = mesh.cell_states[i].ξ

        for j in 1:NC
            mass[j] += ξ[j] * V
        end
    end

    return NTuple{NC, Float64}(mass)
end

"""
    update_mass_tracker!(tracker::MassConservationTracker{NC},
                         J_boundary::Matrix{Float64},
                         A::Float64, dt::Float64) where NC

Update tracker with boundary mass fluxes.

# Arguments
- `tracker`: The mass conservation tracker
- `J_boundary`: Boundary mass fluxes [NC × 2] for top/bottom boundaries [kg/(m²·s)]
- `A`: Cross-sectional area [m²]
- `dt`: Time step [s]

# Notes
- Convention: Positive J = flux into domain (inflow)
- Mass leaving = -J × A × dt (negative J = outflow)
- Typical case: J_top < 0 (gas escaping), J_bottom = 0 (impermeable)
"""
function update_mass_tracker!(tracker::MassConservationTracker{NC},
                              J_boundary::Matrix{Float64},
                              A::Float64, dt::Float64) where NC
    new_efflux = MVector{NC, Float64}(tracker.cumulative_efflux)

    @inbounds for j in 1:NC
        # Sum fluxes from top (index 1) and bottom (index 2) boundaries
        J_top = J_boundary[j, 1]
        J_bottom = J_boundary[j, 2]

        # Mass leaving domain = -J × A × dt (positive J = into domain)
        mass_out = -(J_top + J_bottom) * A * dt
        new_efflux[j] += mass_out
    end

    tracker.cumulative_efflux = NTuple{NC, Float64}(new_efflux)
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

# Returns
NamedTuple with:
- `conserved::Bool`: True if mass is conserved within tolerance
- `current_mass::NTuple{NC, Float64}`: Current mass in domain [kg]
- `expected_mass::NTuple{NC, Float64}`: Expected mass (initial - efflux) [kg]
- `relative_error::NTuple{NC, Float64}`: Relative error per component
- `max_error::Float64`: Maximum relative error across components
"""
function check_mass_conservation(tracker::MassConservationTracker{NC},
                                 mesh::Geometry.Mesh1D{NC}, material;
                                 tol::Float64=1e-10) where NC
    current_mass = compute_domain_mass(mesh, material)

    # Expected mass = initial - cumulative efflux
    expected_mass = ntuple(j -> tracker.initial_mass[j] - tracker.cumulative_efflux[j], Val(NC))

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
        max_error = max_error
    )
end

"""
    mass_conservation_report(tracker::MassConservationTracker{NC},
                             mesh::Mesh1D{NC}, material) -> String

Generate a human-readable mass conservation report.

Useful for post-processing and debugging.
"""
function mass_conservation_report(tracker::MassConservationTracker{NC},
                                  mesh::Geometry.Mesh1D{NC}, material) where NC
    result = check_mass_conservation(tracker, mesh, material)

    lines = String[]
    push!(lines, "=" ^ 60)
    push!(lines, "MASS CONSERVATION REPORT")
    push!(lines, "=" ^ 60)
    push!(lines, @sprintf("Time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 60)

    for j in 1:NC
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

Tracks cumulative energy terms throughout simulation.

For the advection form of convective heat flux, energy conservation is:

    E_current(t) = E_initial + Q_boundary_integrated + Q_reaction_integrated

Where:
- E = ∫ ρcₚT dV  (sensible heat)
- Q_boundary = ∫_∂Ω q·n dA  (boundary heat fluxes: conduction + radiation)
- Q_reaction = ∫ Q_rxn dV  (reaction heat release)

Note: With the advection form, gas enthalpy leaving is already accounted
for in Q_reaction, so no separate boundary enthalpy flux is needed.

# Fields
- `initial_energy::Float64`: Initial sensible heat [J]
- `cumulative_boundary_flux::Float64`: ∫ q_boundary dt × A [J]
- `cumulative_reaction_heat::Float64`: ∫∫ Q_rxn dV dt [J]
- `last_time::Float64`: Time of last update [s]
"""
mutable struct EnergyConservationTracker
    initial_energy::Float64
    cumulative_boundary_flux::Float64
    cumulative_reaction_heat::Float64
    last_time::Float64
end

"""
    EnergyConservationTracker(mesh::Mesh1D, material)

Initialize tracker with current sensible energy in domain.
"""
function EnergyConservationTracker(mesh::Geometry.Mesh1D, material)
    E0 = compute_sensible_energy(mesh, material)
    return EnergyConservationTracker(E0, 0.0, 0.0, 0.0)
end

"""
    compute_sensible_energy(mesh::Mesh1D{NC}, material) -> Float64

Compute total sensible heat in domain: E = Σᵢ (ρcₚ)ᵢ Tᵢ Vᵢ

Only includes active cells (excludes refined-away cells in AMR).

Note: This uses the effective volumetric heat capacity (ρcₚ) computed
from the mixture of solid and gas components at each cell state.
"""
function compute_sensible_energy(mesh::Geometry.Mesh1D{NC}, material) where NC
    E = 0.0

    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end

        state = mesh.cell_states[i]
        T = state.T
        ξ = state.ξ
        V = mesh.cell_volumes[i]

        # Effective ρcₚ at this state
        # This is computed by Properties.effective_heat_capacity
        ρcp = Properties.effective_heat_capacity(ξ, T, material)

        E += ρcp * T * V
    end

    return E
end

"""
    update_energy_tracker!(tracker::EnergyConservationTracker,
                           q_boundary::Float64,
                           Q_rxn_total::Float64,
                           A::Float64, dt::Float64)

Update tracker with boundary fluxes and reaction heat.

# Arguments
- `tracker`: Energy conservation tracker
- `q_boundary`: Net boundary heat flux [W/m²] (positive into domain)
- `Q_rxn_total`: Total reaction heat release ∫ Q_rxn dV [W]
- `A`: Cross-sectional area [m²]
- `dt`: Time step [s]

# Notes
- Boundary flux includes conduction and radiation
- Reaction heat is typically negative for endothermic pyrolysis
- With advection form, gas enthalpy leaving is already in Q_rxn
"""
function update_energy_tracker!(tracker::EnergyConservationTracker,
                                q_boundary::Float64,
                                Q_rxn_total::Float64,
                                A::Float64, dt::Float64)
    # Boundary energy: q × A × dt [J]
    tracker.cumulative_boundary_flux += q_boundary * A * dt

    # Reaction energy: Q × dt [J]
    tracker.cumulative_reaction_heat += Q_rxn_total * dt

    tracker.last_time += dt

    return nothing
end

"""
    check_energy_conservation(tracker::EnergyConservationTracker,
                              mesh::Mesh1D, material;
                              tol::Float64=1e-6) -> NamedTuple

Check energy conservation: E_current = E_initial + Q_boundary + Q_reaction

# Arguments
- `tracker`: Energy conservation tracker
- `mesh`: Current mesh state
- `material`: Material definition
- `tol`: Relative error tolerance (default: 1e-6)

# Returns
NamedTuple with:
- `conserved::Bool`: True if energy is conserved within tolerance
- `current_energy::Float64`: Current sensible heat [J]
- `expected_energy::Float64`: Expected energy from balance [J]
- `initial_energy::Float64`: Initial sensible heat [J]
- `boundary_contribution::Float64`: Cumulative boundary flux [J]
- `reaction_contribution::Float64`: Cumulative reaction heat [J]
- `relative_error::Float64`: Relative error
"""
function check_energy_conservation(tracker::EnergyConservationTracker,
                                   mesh::Geometry.Mesh1D, material;
                                   tol::Float64=1e-6)
    E_current = compute_sensible_energy(mesh, material)
    E_expected = tracker.initial_energy +
                 tracker.cumulative_boundary_flux +
                 tracker.cumulative_reaction_heat

    rel_error = if abs(tracker.initial_energy) > eps()
        abs(E_current - E_expected) / abs(tracker.initial_energy)
    else
        abs(E_current - E_expected)
    end

    return (
        conserved = rel_error < tol,
        current_energy = E_current,
        expected_energy = E_expected,
        initial_energy = tracker.initial_energy,
        boundary_contribution = tracker.cumulative_boundary_flux,
        reaction_contribution = tracker.cumulative_reaction_heat,
        relative_error = rel_error
    )
end

"""
    energy_conservation_report(tracker::EnergyConservationTracker,
                               mesh::Mesh1D, material) -> String

Generate a human-readable energy conservation report.

Useful for post-processing and debugging.
"""
function energy_conservation_report(tracker::EnergyConservationTracker,
                                    mesh::Geometry.Mesh1D, material)
    result = check_energy_conservation(tracker, mesh, material)

    lines = String[]
    push!(lines, "=" ^ 60)
    push!(lines, "ENERGY CONSERVATION REPORT")
    push!(lines, "=" ^ 60)
    push!(lines, @sprintf("Time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 60)
    push!(lines, @sprintf("Initial energy:      %12.6e J", result.initial_energy))
    push!(lines, @sprintf("Boundary flux (Σ):   %12.6e J", result.boundary_contribution))
    push!(lines, @sprintf("Reaction heat (Σ):   %12.6e J", result.reaction_contribution))
    push!(lines, @sprintf("Expected energy:     %12.6e J", result.expected_energy))
    push!(lines, @sprintf("Current energy:      %12.6e J", result.current_energy))
    push!(lines, "-" ^ 60)
    push!(lines, @sprintf("Relative error:      %12.6e", result.relative_error))
    push!(lines, "-" ^ 60)
    status = result.conserved ? "PASSED" : "FAILED"
    push!(lines, @sprintf("Overall: %s", status))
    push!(lines, "=" ^ 60)

    return join(lines, "\n")
end

# ============================================================================
# Per-Reaction Mass Balance Diagnostics
# ============================================================================

"""
    ReactionMassFlowTracker{NC,NR}

Tracks mass flow through each reaction over time.

For each reaction i, tracks:
- Cumulative reactant consumed: Σ∫ rᵢ × νᵢⱼ dt for reactants
- Cumulative product generated: Σ∫ rᵢ × νᵢⱼ dt for products

This allows validation that each reaction conserves mass independently.

# Fields
- `reactant_consumed::Matrix{Float64}`: [NC × NR] cumulative consumption [kg]
- `product_generated::Matrix{Float64}`: [NC × NR] cumulative generation [kg]
- `reaction_counts::Vector{Float64}`: [NR] cumulative ∫ rᵢ dt (extent of reaction) [kg]
- `last_time::Float64`: Time of last update [s]
"""
mutable struct ReactionMassFlowTracker{NC,NR}
    reactant_consumed::Matrix{Float64}   # [NC × NR]
    product_generated::Matrix{Float64}   # [NC × NR]
    reaction_counts::Vector{Float64}     # [NR]
    last_time::Float64
end

"""
    ReactionMassFlowTracker(material::Material{NC,NR}) where {NC,NR}

Initialize reaction mass flow tracker.
"""
function ReactionMassFlowTracker(::Type{<:Material{NC,NR}}) where {NC,NR}
    return ReactionMassFlowTracker{NC,NR}(
        zeros(NC, NR),
        zeros(NC, NR),
        zeros(NR),
        0.0
    )
end

# Convenience constructor that takes a material instance
ReactionMassFlowTracker(material::Material{NC,NR}) where {NC,NR} =
    ReactionMassFlowTracker(typeof(material))

"""
    update_reaction_tracker!(tracker::ReactionMassFlowTracker{NC,NR},
                             rates::Vector{Float64},
                             mesh::Mesh1D{NC},
                             material::Material{NC,NR},
                             dt::Float64) where {NC,NR}

Update reaction mass flow tracker with current rates.

# Arguments
- `tracker`: ReactionMassFlowTracker
- `rates`: Current reaction rates [NR] in kg/(m³·s) (domain-averaged or from a specific cell)
- `mesh`: Mesh (to compute total domain volume)
- `material`: Material definition
- `dt`: Time step [s]

# Notes
This function should be called once per time step with domain-integrated rates.
"""
function update_reaction_tracker!(tracker::ReactionMassFlowTracker{NC,NR},
                                   rates::Vector{Float64},
                                   mesh::Geometry.Mesh1D{NC},
                                   material::Material{NC,NR},
                                   dt::Float64) where {NC,NR}
    # Compute total active volume
    V_total = 0.0
    @inbounds for i in eachindex(mesh.cells)
        if mesh.active_cell_mask[i]
            V_total += mesh.cell_volumes[i]
        end
    end

    @inbounds for i in 1:NR
        rxn = material.reactions[i]
        r = rates[i]
        mass_flow = r * V_total * dt  # kg reacted

        tracker.reaction_counts[i] += mass_flow

        # Track consumption
        if rxn.reactant1 > 0
            tracker.reactant_consumed[rxn.reactant1, i] += rxn.stoich[1] * mass_flow
        end
        if rxn.reactant2 > 0
            tracker.reactant_consumed[rxn.reactant2, i] += rxn.stoich[2] * mass_flow
        end

        # Track generation
        if rxn.product1 > 0
            tracker.product_generated[rxn.product1, i] += rxn.stoich[3] * mass_flow
        end
        if rxn.product2 > 0
            tracker.product_generated[rxn.product2, i] += rxn.stoich[4] * mass_flow
        end
    end

    tracker.last_time += dt
    return nothing
end

"""
    check_reaction_mass_balance(tracker::ReactionMassFlowTracker{NC,NR};
                                 tol::Float64=1e-10) where {NC,NR} -> NamedTuple

Check mass balance for each reaction.

# Arguments
- `tracker`: ReactionMassFlowTracker
- `tol`: Relative tolerance for mass balance error

# Returns
NamedTuple with:
- `balanced::Vector{Bool}`: Whether each reaction is balanced
- `errors::Vector{Float64}`: Relative mass error per reaction
- `all_balanced::Bool`: True if all reactions are balanced
"""
function check_reaction_mass_balance(tracker::ReactionMassFlowTracker{NC,NR};
                                     tol::Float64=1e-10) where {NC,NR}
    balanced = Vector{Bool}(undef, NR)
    errors = Vector{Float64}(undef, NR)

    @inbounds for i in 1:NR
        consumed = sum(tracker.reactant_consumed[:, i])
        generated = sum(tracker.product_generated[:, i])

        if consumed > eps()
            errors[i] = abs(generated - consumed) / consumed
        else
            errors[i] = abs(generated - consumed)
        end

        balanced[i] = errors[i] < tol
    end

    return (balanced=balanced, errors=errors, all_balanced=all(balanced))
end

"""
    reaction_mass_flow_report(tracker::ReactionMassFlowTracker{NC,NR},
                              material::Material{NC,NR}) where {NC,NR} -> String

Generate detailed report of mass flow through reactions.

# Example
```julia
tracker = ReactionMassFlowTracker(material)
# ... run simulation with update_reaction_tracker! calls ...
println(reaction_mass_flow_report(tracker, material))
```
"""
function reaction_mass_flow_report(tracker::ReactionMassFlowTracker{NC,NR},
                                   material::Material{NC,NR}) where {NC,NR}
    result = check_reaction_mass_balance(tracker)

    lines = String[]
    push!(lines, "=" ^ 70)
    push!(lines, "REACTION MASS FLOW REPORT")
    push!(lines, "=" ^ 70)
    push!(lines, @sprintf("Simulation time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 70)

    for i in 1:NR
        status = result.balanced[i] ? "✓" : "✗"
        push!(lines, @sprintf("\nReaction %d %s (error: %.2e)", i, status, result.errors[i]))

        consumed = sum(tracker.reactant_consumed[:, i])
        generated = sum(tracker.product_generated[:, i])

        push!(lines, @sprintf("  Total consumed:  %.6e kg", consumed))
        push!(lines, @sprintf("  Total generated: %.6e kg", generated))

        # Per-component breakdown
        push!(lines, "  By component:")
        for j in 1:NC
            c = tracker.reactant_consumed[j, i]
            g = tracker.product_generated[j, i]
            if c > 0 || g > 0
                comp_name = string(material.components[j].name)
                push!(lines, @sprintf("    %-12s: consumed=%.4e, generated=%.4e", comp_name, c, g))
            end
        end
    end

    push!(lines, "\n" * "-" ^ 70)
    overall = result.all_balanced ? "PASSED" : "FAILED"
    push!(lines, @sprintf("Overall: %s", overall))
    push!(lines, "=" ^ 70)

    return join(lines, "\n")
end
