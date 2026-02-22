"""
Time step control for explicit and semi-implicit integrators.

Reference: Design Document Section 8.3

This module provides CFL-based timestep estimation for stability control.
The timestep is determined by the minimum of diffusive and advective limits.
"""

using ..Geometry: cell_size, face_position

# ============================================================================
# Timestep Computation
# ============================================================================

"""
    compute_timestep(mesh, prop_cache::PropertyCache;
                     cfl::Float64=0.5, safety::Float64=0.9) -> Float64

Compute a stable timestep based on CFL conditions.

The timestep is the minimum of:
- Diffusive limit: dt_diff = Δx² / (2α) where α = k/(ρcₚ)
- Advective limit: dt_adv = Δx / |v| (if applicable)

# Arguments
- `mesh`: Mesh with geometry information
- `prop_cache`: Property cache with evaluated properties
- `cfl`: CFL number (default 0.5)
- `safety`: Additional safety factor (default 0.9)

# Returns
Stable timestep Δt [s]
"""
function compute_timestep(mesh, prop_cache::Properties.PropertyCache;
                          cfl::Float64=0.5, safety::Float64=0.9)
    n_cells_val = length(mesh.cells)

    if n_cells_val == 0
        return 1.0  # No cells, return arbitrary timestep
    end

    dt_min = Inf

    @inbounds for i in 1:n_cells_val
        # Get cell size
        Δx = cell_size(mesh, i)

        # Get thermal properties
        k = prop_cache.k[i]
        ρcp = prop_cache.ρcp[i]

        # Thermal diffusivity
        if ρcp > 0.0 && k > 0.0
            α = k / ρcp

            # Diffusive stability limit: dt < Δx²/(2α)
            dt_diff = Δx^2 / (2.0 * α)
            dt_min = min(dt_min, dt_diff)
        end
    end

    # Apply CFL and safety factors
    dt = cfl * safety * dt_min

    # Guard against zero or negative timestep
    if dt <= 0.0 || !isfinite(dt)
        dt = 1e-10
    end

    return dt
end

"""
    compute_timestep(mesh, material, T::Float64=300.0;
                     cfl::Float64=0.5, safety::Float64=0.9) -> Float64

Compute timestep using constant temperature for property evaluation.

Useful for initial timestep estimation before solving.

# Arguments
- `mesh`: Mesh with geometry
- `material`: Material for property evaluation
- `T`: Temperature for property evaluation [K] (default 300)
- `cfl`: CFL number
- `safety`: Safety factor
"""
function compute_timestep(mesh, material, T::Float64=300.0;
                          cfl::Float64=0.5, safety::Float64=0.9)
    n_cells_val = length(mesh.cells)

    if n_cells_val == 0
        return 1.0
    end

    dt_min = Inf

    @inbounds for i in 1:n_cells_val
        Δx = cell_size(mesh, i)
        state = mesh.cell_states[i]

        # Compute effective properties at current state
        k = Properties.effective_conductivity(state.ξ, T, material)
        ρcp = Properties.effective_heat_capacity(state.ξ, T, material)

        if ρcp > 0.0 && k > 0.0
            α = k / ρcp
            dt_diff = Δx^2 / (2.0 * α)
            dt_min = min(dt_min, dt_diff)
        end
    end

    dt = cfl * safety * dt_min

    if dt <= 0.0 || !isfinite(dt)
        dt = 1e-10
    end

    return dt
end

# ============================================================================
# Reaction-Limited Timestep
# ============================================================================

"""
    compute_reaction_timestep(mesh, material;
                              max_fraction::Float64=0.1) -> Float64

Estimate timestep limit from reaction rates.

Ensures that no more than `max_fraction` of any component is consumed
in a single timestep: dt < max_fraction * ξ / |dξ/dt|

# Arguments
- `mesh`: Mesh with current state
- `material`: Material with reactions
- `max_fraction`: Maximum fraction change per step (default 0.1)

# Returns
Reaction-limited timestep [s]
"""
function compute_reaction_timestep(mesh, material;
                                   max_fraction::Float64=0.1)
    n_cells_val = length(mesh.cells)
    NR = n_reactions(material)

    if n_cells_val == 0 || NR == 0
        return Inf
    end

    dt_min = Inf
    rates = zeros(NR)
    heats = zeros(NR)

    @inbounds for i in 1:n_cells_val
        state = mesh.cell_states[i]
        T = state.T
        ξ = state.ξ

        # Evaluate reaction rates
        Physics.evaluate_reactions!(rates, heats, ξ, T, material)

        # Check each reactant consumption rate
        for (r, rxn) in enumerate(material.reactions)
            rate = rates[r]
            if rate > 0.0
                # Check reactant 1
                if rxn.reactant1 > 0
                    ξ_r = ξ[rxn.reactant1]
                    if ξ_r > 0.0
                        # Rate of consumption (negative)
                        dξ_dt = -rxn.stoich[1] * rate
                        if dξ_dt < 0.0  # Consumption
                            dt_rxn = max_fraction * ξ_r / abs(dξ_dt)
                            dt_min = min(dt_min, dt_rxn)
                        end
                    end
                end

                # Check reactant 2
                if rxn.reactant2 > 0
                    ξ_r = ξ[rxn.reactant2]
                    if ξ_r > 0.0
                        dξ_dt = -rxn.stoich[2] * rate
                        if dξ_dt < 0.0
                            dt_rxn = max_fraction * ξ_r / abs(dξ_dt)
                            dt_min = min(dt_min, dt_rxn)
                        end
                    end
                end
            end
        end
    end

    return dt_min
end

# ============================================================================
# Combined Timestep
# ============================================================================

"""
    compute_stable_timestep(mesh, prop_cache::PropertyCache, material;
                            cfl::Float64=0.5,
                            safety::Float64=0.9,
                            max_reaction_fraction::Float64=0.1,
                            dt_max::Float64=Inf,
                            dt_min::Float64=1e-15) -> Float64

Compute a stable timestep considering all constraints.

Takes the minimum of:
- Diffusive CFL limit
- Reaction rate limit
- User-specified maximum

# Arguments
- `mesh`: Mesh with state
- `prop_cache`: Property cache
- `material`: Material definition
- `cfl`: CFL number for diffusion
- `safety`: Safety factor
- `max_reaction_fraction`: Max species change fraction
- `dt_max`: Maximum allowed timestep
- `dt_min`: Minimum allowed timestep

# Returns
Stable timestep [s]
"""
function compute_stable_timestep(mesh, prop_cache::Properties.PropertyCache, material;
                                 cfl::Float64=0.5,
                                 safety::Float64=0.9,
                                 max_reaction_fraction::Float64=0.1,
                                 dt_max::Float64=Inf,
                                 dt_min::Float64=1e-15)
    # Diffusive limit
    dt_diff = compute_timestep(mesh, prop_cache; cfl=cfl, safety=safety)

    # Reaction limit
    dt_rxn = compute_reaction_timestep(mesh, material; max_fraction=max_reaction_fraction)

    # Combined limit
    dt = min(dt_diff, dt_rxn, dt_max)

    # Apply floor
    dt = max(dt, dt_min)

    return dt
end

# ============================================================================
# Timestep Controller for Adaptive Stepping
# ============================================================================

"""
    TimestepController

Adaptive timestep controller with error-based adjustment.

# Fields
- `dt::Float64`: Current timestep
- `dt_min::Float64`: Minimum allowed timestep
- `dt_max::Float64`: Maximum allowed timestep
- `growth_factor::Float64`: Maximum growth per step
- `shrink_factor::Float64`: Minimum shrink per step
- `target_error::Float64`: Target local error
"""
mutable struct TimestepController
    dt::Float64
    dt_min::Float64
    dt_max::Float64
    growth_factor::Float64
    shrink_factor::Float64
    target_error::Float64
end

"""
    TimestepController(; dt_initial=1e-6, dt_min=1e-15, dt_max=1.0,
                        growth_factor=1.5, shrink_factor=0.5,
                        target_error=1e-4)

Create a timestep controller with default parameters.
"""
function TimestepController(; dt_initial::Float64=1e-6,
                              dt_min::Float64=1e-15,
                              dt_max::Float64=1.0,
                              growth_factor::Float64=1.5,
                              shrink_factor::Float64=0.5,
                              target_error::Float64=1e-4)
    TimestepController(dt_initial, dt_min, dt_max,
                       growth_factor, shrink_factor, target_error)
end

"""
    adjust_timestep!(controller::TimestepController, error::Float64) -> Float64

Adjust timestep based on error estimate.

# Arguments
- `controller`: Timestep controller
- `error`: Estimated local error

# Returns
New timestep
"""
function adjust_timestep!(controller::TimestepController, error::Float64)
    if error <= 0.0
        # No error info, grow cautiously
        controller.dt = min(controller.dt * controller.growth_factor, controller.dt_max)
    elseif error > controller.target_error
        # Error too large, shrink
        factor = (controller.target_error / error)^0.5
        factor = max(factor, controller.shrink_factor)
        controller.dt = max(controller.dt * factor, controller.dt_min)
    else
        # Error acceptable, try to grow
        factor = (controller.target_error / error)^0.5
        factor = min(factor, controller.growth_factor)
        controller.dt = min(controller.dt * factor, controller.dt_max)
    end

    return controller.dt
end

"""
    reset_timestep!(controller::TimestepController, dt::Float64)

Reset controller to a specific timestep.
"""
function reset_timestep!(controller::TimestepController, dt::Float64)
    controller.dt = clamp(dt, controller.dt_min, controller.dt_max)
    return controller.dt
end
