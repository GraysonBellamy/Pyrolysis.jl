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

# ALE gas-species note
The gas mesh-motion correction is the **conservative flux form** (upwinded
face sweep fluxes `[ξ_g w A]` plus the `ξ_g V̇/V` volume-rate term,
`compute_ale_species_advection!`), so the discrete gas ledger telescopes to
boundary fluxes exactly — the same bookkeeping-level closure as condensed
components. The `expected_mass` formula here counts only the *transport*
efflux (`J_gas` at the boundary faces): under ALE with an **open** top/bottom
mass BC, the moving boundary additionally sweeps gas in/out at
`ξ_face·w_face·A` (a regressing surface releases the pore gas it uncovers),
which this tracker does not yet integrate — expect a corresponding offset in
open-boundary ALE runs. Fully closed (impermeable) columns conserve each gas
component to integrator tolerance (V-CPL-03 P4).

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

Tracks **two complementary energy ledgers** for the Formulation A
pyrolysis residual:

## 1. Discrete-form ledger (closes to integrator tolerance)

Built from the terms the discretized RHS *actually* assembles:

    ΔE_matrix(t) = ∫ q_cond_boundary · A · dt
                 + ∫ Q_rxn dV · dt
                 + ∫ S_conv dV · dt
                 + ∫ S_gen dV · dt          (zero unless `:with_generation_sink`)

with `E_matrix = Σ_i ρcp_eff,i · T_i · V_i` using the **gas-excluded**
`ρcp_eff` (quasi-steady-gas approximation; NOT the ThermaKin convention, which
keeps gas storage — see `Properties.effective_heat_capacity`). This residual
closes to roughly integrator tolerance × |E(0)| by construction; failures here
indicate bugs in flux/source assembly or callbacks. NOTE: because it is built
from the terms the RHS *actually* assembles, this ledger closes even if the RHS
itself double-counts a mechanism — it audits discrete bookkeeping, not the
physical model. `cumulative_transpiration` is zero under the default
(transpiration-off) convention.

## 2. Physical-form ledger (auditable cost of the FDS-convention approximation)

Built from total sensible energy including gas storage. The physical-form
residual quantifies the cost of the quasi-steady-gas approximation:

    E_total = Σ_i ρcp_total,i · T_i · V_i     (matrix + gas)
    physical_residual(t) = E_total(t) − E_total(0)
                          − cumulative_boundary_applied
                          − cumulative_reaction_heat
                          + cumulative_transpiration

For canonical pyrolysis cases (gas residence time ≪ thermal time scale)
this should be < 1% × |E(0)|; if it is larger, the regime stresses the
quasi-steady-gas approximation and the user should reconsider.

# Fields
- `initial_energy_matrix::Float64`     — `Σ ρcp_eff · T · V` at t=0 (matrix only).
- `initial_energy_total::Float64`      — `Σ ρcp_total · T · V` at t=0 (matrix + gas).
- `cumulative_boundary_conduction::Float64` — `∫ q_cond_boundary · A · dt`.
- `cumulative_boundary_applied::Float64`    — `∫ q_BC_applied · A · dt`.
- `cumulative_reaction_heat::Float64`       — `∫ Q_rxn dV · dt`.
- `cumulative_advection_sink::Float64`      — `∫ ∫ S_conv dV · dt`.
- `cumulative_generation_sink::Float64`    — `∫ ∫ S_gen dV · dt`.
- `cumulative_transpiration::Float64`       — `∫ Σ_g ṁ_g · Δh_g(T_s, T_∞) · A · dt`.
- `last_time::Float64`                      — Time of last update [s].
"""
mutable struct EnergyConservationTracker
    initial_energy_matrix::Float64
    initial_energy_total::Float64
    cumulative_boundary_conduction::Float64
    cumulative_boundary_applied::Float64
    cumulative_reaction_heat::Float64
    cumulative_advection_sink::Float64
    cumulative_generation_sink::Float64
    cumulative_transpiration::Float64
    last_time::Float64
end

"""
    EnergyConservationTracker(mesh::Mesh1D, material)

Initialize tracker with current sensible energies (matrix-only and
total-including-gas).
"""
function EnergyConservationTracker(mesh::Geometry.Mesh1D, material)
    E0_matrix = compute_sensible_energy(mesh, material)
    E0_total = compute_total_sensible_energy(mesh, material)
    return EnergyConservationTracker(E0_matrix, E0_total, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    compute_sensible_energy(mesh::Mesh1D{NC}, material) -> Float64

Matrix-only sensible heat: `E_matrix = Σ_i ρcp_eff,i · T_i · V_i`. Gas
components are excluded from `ρcp_eff` per Formulation A. Only active cells
(AMR-aware).
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
    compute_total_sensible_energy(mesh::Mesh1D{NC}, material) -> Float64

Total sensible heat (matrix + gas storage):
`E_total = Σ_i ρcp_total,i · T_i · V_i`, where `ρcp_total` includes all
components per `Properties.total_heat_capacity`. Used by the physical-form
energy residual to quantify the cost of the quasi-steady-gas approximation
that the matrix-only `ρcp_eff` makes.
"""
function compute_total_sensible_energy(mesh::Geometry.Mesh1D{NC}, material) where {NC}
    E = 0.0
    @inbounds for i in eachindex(mesh.cells)
        if !mesh.active_cell_mask[i]
            continue
        end

        state = mesh.cell_states[i]
        T = state.T
        ξ = state.ξ
        V = mesh.cell_volumes[i]

        ρcp = Properties.total_heat_capacity(ξ, T, material)
        E += ρcp * T * V
    end
    return E
end

"""
    update_energy_tracker!(tracker; q_cond_boundary, q_bc_applied,
                           Q_rxn_total, S_conv_total, S_gen_total,
                           q_transpiration, A, dt)

Update tracker with one substep's contributions.

All `*_total` arguments are **volume-integrals over the cell-source field**
(`Σ_i field_i · V_i`) computed from the active RHS state. `q_*_boundary`
are face-flux values [W/m²] (multiplied by A and dt internally).

# Keyword arguments
- `q_cond_boundary::Float64` — net conducted boundary flux [W/m²], positive into slab.
- `q_bc_applied::Float64`    — applied external BC flux [W/m²], positive into slab.
- `Q_rxn_total::Float64`     — `∫ Q_rxn dV` [W].
- `S_conv_total::Float64`    — `∫ S_conv dV` [W] (typically negative).
- `S_gen_total::Float64`   — `∫ S_gen dV` [W] (zero unless `:with_generation_sink`).
- `q_transpiration::Float64` — `Σ_g ṁ_g · Δh_g(T_s, T_∞)` summed over boundary faces [W/m²].
- `A::Float64`               — Cross-sectional area at this substep [m²].
- `dt::Float64`              — Substep size [s].
"""
function update_energy_tracker!(
    tracker::EnergyConservationTracker;
    q_cond_boundary::Float64,
    q_bc_applied::Float64,
    Q_rxn_total::Float64,
    S_conv_total::Float64,
    S_gen_total::Float64,
    q_transpiration::Float64,
    A::Float64,
    dt::Float64,
)
    tracker.cumulative_boundary_conduction += q_cond_boundary * A * dt
    tracker.cumulative_boundary_applied += q_bc_applied * A * dt
    tracker.cumulative_reaction_heat += Q_rxn_total * dt
    tracker.cumulative_advection_sink += S_conv_total * dt
    tracker.cumulative_generation_sink += S_gen_total * dt
    tracker.cumulative_transpiration += q_transpiration * A * dt
    tracker.last_time += dt
    return nothing
end

"""
    check_energy_conservation(tracker, mesh, material; tol=1e-6) -> NamedTuple

Compute both the discrete-form residual (closes to integrator tolerance by
construction) and the physical-form residual (quantifies the quasi-steady-gas
approximation cost).

# Returns
NamedTuple with:
- `discrete_residual::Float64`        — `E_matrix(t) − E_matrix(0) − Σ contributions`. Should be `≈ 0` to integrator tolerance.
- `discrete_relative::Float64`        — relative magnitude vs `|E_matrix(0)|`.
- `physical_residual::Float64`        — `E_total(t) − E_total(0) − applied + transpiration − reaction`.
- `physical_relative::Float64`        — relative magnitude vs `|E_total(0)|`.
- `current_energy_matrix::Float64`    — current matrix sensible heat.
- `current_energy_total::Float64`     — current total sensible heat (matrix + gas).
- `gas_storage_gap::Float64`          — `(E_total − E_matrix)`, the energy held in gas storage that the discrete equation does not evolve.
- `discrete_conserved::Bool`          — `discrete_relative < tol`.
- `physical_conserved::Bool`          — `physical_relative < tol`.
"""
function check_energy_conservation(
    tracker::EnergyConservationTracker,
    mesh::Geometry.Mesh1D,
    material;
    tol::Float64 = 1e-6,
)
    E_matrix_now = compute_sensible_energy(mesh, material)
    E_total_now = compute_total_sensible_energy(mesh, material)

    # Discrete-form expectation: only terms the RHS actually assembles
    # (matrix-only ρcp; S_conv and S_gen already integrated).
    expected_matrix =
        tracker.initial_energy_matrix +
        tracker.cumulative_boundary_conduction +
        tracker.cumulative_reaction_heat +
        tracker.cumulative_advection_sink +
        tracker.cumulative_generation_sink
    discrete_residual = E_matrix_now - expected_matrix

    # Physical-form expectation: full energy ledger including gas storage,
    # using applied BC heat (not just conducted) and the transpiration outflow
    # as the explicit energy "leaving via gas".
    expected_total =
        tracker.initial_energy_total +
        tracker.cumulative_boundary_applied +
        tracker.cumulative_reaction_heat - tracker.cumulative_transpiration
    physical_residual = E_total_now - expected_total

    discrete_relative = if abs(tracker.initial_energy_matrix) > eps()
        abs(discrete_residual) / abs(tracker.initial_energy_matrix)
    else
        abs(discrete_residual)
    end

    physical_relative = if abs(tracker.initial_energy_total) > eps()
        abs(physical_residual) / abs(tracker.initial_energy_total)
    else
        abs(physical_residual)
    end

    return (
        discrete_residual = discrete_residual,
        discrete_relative = discrete_relative,
        discrete_conserved = discrete_relative < tol,
        physical_residual = physical_residual,
        physical_relative = physical_relative,
        physical_conserved = physical_relative < tol,
        current_energy_matrix = E_matrix_now,
        current_energy_total = E_total_now,
        gas_storage_gap = E_total_now - E_matrix_now,
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
    push!(lines, "ENERGY CONSERVATION REPORT (Formulation A)")
    push!(lines, "=" ^ 70)
    push!(lines, @sprintf("Time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 70)
    push!(
        lines,
        "DISCRETE-FORM LEDGER (matrix-only ρcp_eff; closes to integrator tolerance)",
    )
    push!(
        lines,
        @sprintf("  E_matrix(0):                %12.6e J", tracker.initial_energy_matrix)
    )
    push!(
        lines,
        @sprintf("  E_matrix(t):                %12.6e J", result.current_energy_matrix)
    )
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
            "  ∫∫ S_gen dV · dt:         %12.6e J",
            tracker.cumulative_generation_sink
        )
    )
    push!(
        lines,
        @sprintf(
            "  Discrete residual:          %12.6e J  (rel. %.2e) %s",
            result.discrete_residual,
            result.discrete_relative,
            result.discrete_conserved ? "✓" : "✗"
        )
    )
    push!(lines, "-" ^ 70)
    push!(
        lines,
        "PHYSICAL-FORM LEDGER (matrix + gas storage; quantifies approximation cost)",
    )
    push!(
        lines,
        @sprintf("  E_total(0):                 %12.6e J", tracker.initial_energy_total)
    )
    push!(
        lines,
        @sprintf("  E_total(t):                 %12.6e J", result.current_energy_total)
    )
    push!(lines, @sprintf("  Gas storage gap:            %12.6e J", result.gas_storage_gap))
    push!(
        lines,
        @sprintf(
            "  ∫ q_BC_applied · A · dt:    %12.6e J",
            tracker.cumulative_boundary_applied
        )
    )
    push!(
        lines,
        @sprintf(
            "  ∫ Σ ṁ·Δh · A · dt (transp): %12.6e J",
            tracker.cumulative_transpiration
        )
    )
    push!(
        lines,
        @sprintf(
            "  Physical residual:          %12.6e J  (rel. %.2e) %s",
            result.physical_residual,
            result.physical_relative,
            result.physical_conserved ? "✓" : "✗"
        )
    )
    push!(lines, "=" ^ 70)
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
    return ReactionMassFlowTracker{NC,NR}(zeros(NC, NR), zeros(NC, NR), zeros(NR), 0.0)
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
function update_reaction_tracker!(
    tracker::ReactionMassFlowTracker{NC,NR},
    rates::Vector{Float64},
    mesh::Geometry.Mesh1D{NC},
    material::Material{NC,NR},
    dt::Float64,
) where {NC,NR}
    # Compute total active volume
    V_total = 0.0
    @inbounds for i in eachindex(mesh.cells)
        if mesh.active_cell_mask[i]
            V_total += mesh.cell_volumes[i]
        end
    end

    @inbounds for i = 1:NR
        rxn = material.reactions[i]
        r = rates[i]
        mass_flow = r * V_total * dt  # kg reacted

        tracker.reaction_counts[i] += mass_flow

        # Track consumption: convention — first reactant consumes mass 1·r per unit extent
        if length(rxn.reactants) >= 1
            r1_idx, _ = rxn.reactants[1]
            tracker.reactant_consumed[r1_idx, i] += mass_flow
        end

        # Track generation: each product contributes its mass yield
        for (p_idx, ν) in rxn.products
            tracker.product_generated[p_idx, i] += ν * mass_flow
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
function check_reaction_mass_balance(
    tracker::ReactionMassFlowTracker{NC,NR};
    tol::Float64 = 1e-10,
) where {NC,NR}
    balanced = Vector{Bool}(undef, NR)
    errors = Vector{Float64}(undef, NR)

    @inbounds for i = 1:NR
        consumed = sum(tracker.reactant_consumed[:, i])
        generated = sum(tracker.product_generated[:, i])

        if consumed > eps()
            errors[i] = abs(generated - consumed) / consumed
        else
            errors[i] = abs(generated - consumed)
        end

        balanced[i] = errors[i] < tol
    end

    return (balanced = balanced, errors = errors, all_balanced = all(balanced))
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
function reaction_mass_flow_report(
    tracker::ReactionMassFlowTracker{NC,NR},
    material::Material{NC,NR},
) where {NC,NR}
    result = check_reaction_mass_balance(tracker)

    lines = String[]
    push!(lines, "=" ^ 70)
    push!(lines, "REACTION MASS FLOW REPORT")
    push!(lines, "=" ^ 70)
    push!(lines, @sprintf("Simulation time: %.4f s", tracker.last_time))
    push!(lines, "-" ^ 70)

    for i = 1:NR
        status = result.balanced[i] ? "✓" : "✗"
        push!(
            lines,
            @sprintf("\nReaction %d %s (error: %.2e)", i, status, result.errors[i])
        )

        consumed = sum(tracker.reactant_consumed[:, i])
        generated = sum(tracker.product_generated[:, i])

        push!(lines, @sprintf("  Total consumed:  %.6e kg", consumed))
        push!(lines, @sprintf("  Total generated: %.6e kg", generated))

        # Per-component breakdown
        push!(lines, "  By component:")
        for j = 1:NC
            c = tracker.reactant_consumed[j, i]
            g = tracker.product_generated[j, i]
            if c > 0 || g > 0
                comp_name = string(material.components[j].name)
                push!(
                    lines,
                    @sprintf("    %-12s: consumed=%.4e, generated=%.4e", comp_name, c, g)
                )
            end
        end
    end

    push!(lines, "\n" * "-" ^ 70)
    overall = result.all_balanced ? "PASSED" : "FAILED"
    push!(lines, @sprintf("Overall: %s", overall))
    push!(lines, "=" ^ 70)

    return join(lines, "\n")
end
