"""
Finite volume operators.

This module provides:
- Divergence operator for 1D finite volumes
- RHS assembly for the complete system of equations
"""

using ..Materials: R_GAS, STEFAN_BOLTZMANN, n_components, n_reactions
using ..BoundaryConditions: DirichletThermalBC
using ..Materials: is_gas
using ..Materials: evaluate_param
using ..Geometry: cell_size
using ..Properties: PropertyCache, update_properties!, effective_permeability, porosity
using ..Physics:
    compute_conductive_fluxes!,
    compute_radiation_source!,
    compute_gas_fluxes!,
    evaluate_reactions!,
    species_source_terms!,
    reaction_heat_source

# ============================================================================
# Divergence Operator
# ============================================================================

"""
	divergence_1d!(div::Vector{Float64}, fluxes::Vector{Float64}, mesh::Mesh1D)

Compute the divergence of a face-centered flux field.

For cell i with faces f_left and f_right:
	div[i] = (1/V_i) * (F[f_right] - F[f_left]) * A

Sign convention:
- Positive flux = flow in +z direction (from left to right / top to bottom)
- For a cell, right face flux is outflow (+), left face flux is inflow (-)

# Arguments
- `div`: Output vector of divergence values [units of flux / length]
- `fluxes`: Face-centered flux values [flux units]
- `mesh`: 1D mesh with geometry
"""
function divergence_1d!(div::AbstractVector, fluxes::AbstractVector, state, mesh)
    n_cells_val = mesh.n_active_cells

    @inbounds @simd for i = 1:n_cells_val
        cell = mesh.cells[i]
        f_left, f_right = cell.face_ids

        V_i = state.volume[i]
        A = state.A_face[f_left]

        div[i] = (fluxes[f_right] - fluxes[f_left]) * A / V_i
    end

    return nothing
end

"""
	divergence_species_1d!(div::Matrix{Float64}, J::Matrix{Float64}, mesh::Mesh1D)

Compute divergence for species fluxes (multiple components).

# Arguments  
- `div`: Output matrix [NC × n_cells] of divergence values
- `J`: Face-centered species fluxes [NC × n_faces]
- `mesh`: 1D mesh
"""
function divergence_species_1d!(div::AbstractMatrix, J::AbstractMatrix, state, mesh)
    NC = size(div, 1)
    n_cells_val = mesh.n_active_cells

    @inbounds for i = 1:n_cells_val
        cell = mesh.cells[i]
        f_left, f_right = cell.face_ids

        V_i = state.volume[i]
        A = state.A_face[f_left]

        @simd for j = 1:NC
            div[j, i] = (J[j, f_right] - J[j, f_left]) * A / V_i
        end
    end

    return nothing
end

# ============================================================================
# RHS Assembly for Fixed Mesh (Eulerian)
# ============================================================================

"""
	RHSCache{Tu}

Pre-allocated workspace for RHS computation. Parameterized on the live
scalar type `Tu` so the AD workspace can hold `ForwardDiff.Dual`-typed
flux/source values without re-allocation.
"""
mutable struct RHSCache{Tu}
    # Work arrays for fluxes
    q_cond::Vector{Tu}      # Conductive heat fluxes at faces
    q_total::Vector{Tu}     # Total energy fluxes at faces
    S_conv::Vector{Tu}      # Gas convective volumetric source per cell [W/m³]
    J_gas::Matrix{Tu}       # Gas species fluxes [NC × n_faces]

    # Work arrays for sources
    Q_rxn::Vector{Tu}       # Reaction heat source per cell
    Q_rad::Vector{Tu}       # Radiation source per cell
    dξ_rxn::Matrix{Tu}      # Species source from reactions [NC × n_cells]

    # Per-cell reaction workspace
    rates::Vector{Tu}       # Reaction rates (NR)
    heats::Vector{Tu}       # Heats of reaction (NR)
    species_src::Vector{Tu} # Species sources (NC)

    # Divergence workspace
    div_q::Vector{Tu}       # Divergence of heat flux
    div_J::Matrix{Tu}       # Divergence of species flux [NC × n_cells]
end

"""
	RHSCache{Tu}(n_cells::Int, n_faces::Int, NC::Int, NR::Int) where {Tu}
	RHSCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)         # Tu = Float64

Create an RHSCache with pre-allocated arrays.
"""
function RHSCache{Tu}(n_cells::Int, n_faces::Int, NC::Int, NR::Int) where {Tu}
    RHSCache{Tu}(
        zeros(Tu, n_faces),           # q_cond
        zeros(Tu, n_faces),           # q_total
        zeros(Tu, n_cells),           # S_conv
        zeros(Tu, NC, n_faces),       # J_gas
        zeros(Tu, n_cells),           # Q_rxn
        zeros(Tu, n_cells),           # Q_rad
        zeros(Tu, NC, n_cells),       # dξ_rxn
        zeros(Tu, NR),                # rates
        zeros(Tu, NR),                # heats
        zeros(Tu, NC),                # species_src
        zeros(Tu, n_cells),           # div_q
        zeros(Tu, NC, n_cells),       # div_J
    )
end

RHSCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int) =
    RHSCache{Float64}(n_cells, n_faces, NC, NR)

"""
	RHSCache(mesh, material)

Create a Float64-typed RHSCache sized for the given mesh and material.
"""
function RHSCache(mesh, material)
    n_cells_val = length(mesh.cells)
    n_faces_val = length(mesh.faces)
    NC = n_components(material)
    NR = n_reactions(material)
    return RHSCache(n_cells_val, n_faces_val, NC, NR)
end

"""
	compute_rhs!(dT::Vector{Float64}, dξ::Matrix{Float64},
				 mesh::Mesh1D{NC}, prop_cache::PropertyCache, rhs_cache::RHSCache,
				 material::Material{NC}, bc_set::BoundaryConditionSet, t,
				 top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int};
				 I_radiation::Float64=0.0) where NC

Compute the right-hand side of the governing equations using type-stable BCs.

This version uses pre-cached face lists and static dispatch for maximum
performance. All boundary condition handling is specialized at compile time.

# Arguments
- `dT`: Output vector of dT/dt [K/s]
- `dξ`: Output matrix of dξⱼ/dt [kg/(m³·s)], size [NC × n_cells]
- `mesh`: Mesh with current state (T, ξ in cell_states)
- `prop_cache`: Property cache (must be updated before calling)
- `rhs_cache`: RHS workspace cache
- `material`: Material definition
- `bc_set`: Type-stable BoundaryConditionSet with all boundary conditions
- `t`: Current time [s]
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `I_radiation`: Surface radiation intensity [W/m²] (default 0)

# Equations solved
Energy: ρcₚ dT/dt = -∇·q_cond + S_conv + Q_rad + Q_rxn
		(S_conv = -Σ_g cp_g·J̄_g·∂T/∂x)
Species: dξⱼ/dt = -∇·Jⱼ + Sⱼ (for gases with transport)
		 dξⱼ/dt = Sⱼ          (for solids, no transport)

The `NC` type parameter enables the compiler to specialize for the specific
number of components, eliminating type instability and allocations in the
inner loops.
"""
function compute_rhs!(
    dT::AbstractVector,
    dξ::AbstractMatrix,
    state,
    mesh::Mesh1D{NC},
    prop_cache::PropertyCache,
    rhs_cache::RHSCache,
    material::Material{NC},
    bc_set::BoundaryConditionSet,
    t,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int};
    I_radiation::Real = 0.0,
    radiation_model::RadiationModel = NO_RADIATION,
    radiation_kirchhoff::Bool = true,
) where {NC}

    n_cells_val = mesh.n_active_cells

    # ========================================================================
    # Step 1: Compute conductive heat fluxes at interior faces
    # ========================================================================
    compute_conductive_fluxes!(rhs_cache.q_cond, state, mesh, prop_cache)

    # Copy to total flux (will add boundary contributions)
    copyto!(rhs_cache.q_total, rhs_cache.q_cond)

    # ========================================================================
    # Step 2: Compute reaction rates and source terms for each cell
    # ========================================================================
    @inbounds for i = 1:n_cells_val
        T = state.T[i]
        ξ = ntuple(j -> state.ξ[j, i], Val(NC))

        # Evaluate reaction rates and heats
        evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, ξ, T, material)

        # Compute reaction heat source
        rhs_cache.Q_rxn[i] = reaction_heat_source(rhs_cache.rates, rhs_cache.heats)

        # Compute species source terms
        species_source_terms!(rhs_cache.species_src, rhs_cache.rates, material)

        # Copy to output matrix
        for j = 1:NC
            rhs_cache.dξ_rxn[j, i] = rhs_cache.species_src[j]
        end
    end

    # ========================================================================
    # Step 3: Compute gas transport fluxes (interior + boundary)
    # Must run BEFORE thermal BCs: the blowing models of the convective BCs
    # read the outward gas flux from rhs_cache.J_gas at the boundary face.
    # ========================================================================
    compute_gas_fluxes!(rhs_cache.J_gas, state, mesh, prop_cache, material)

    apply_mass_bc_flux!(
        rhs_cache.J_gas,
        state,
        mesh,
        prop_cache,
        bc_set,
        material,
        t,
        top_face_ids,
        bottom_face_ids,
    )

    # ========================================================================
    # Step 4: Apply thermal boundary conditions (static dispatch)
    # ========================================================================
    apply_thermal_bc_flux!(
        rhs_cache.q_total,
        state,
        mesh,
        bc_set,
        prop_cache,
        material,
        rhs_cache.J_gas,
        t,
        top_face_ids,
        bottom_face_ids,
        radiation_model,
    )

    # ========================================================================
    # Step 5: Compute radiation source (model-dependent)
    # ========================================================================
    compute_radiation_source!(
        rhs_cache.Q_rad,
        state,
        mesh,
        prop_cache,
        I_radiation,
        radiation_model,
        radiation_kirchhoff,
    )

    # ========================================================================
    # Step 6: Cell-centered gas convective source
    # ------------------------------------------------------------------------
    # S_conv = -Σ_g cp_g(T_i) · J̄_g,i · ∂T/∂x|_i
    # Needs q_total: the boundary cells' gradient blends in the boundary-face
    # slope −q/k at the T_s solved in Step 4.
    # ========================================================================
    Physics.compute_gas_convective_source!(
        rhs_cache.S_conv,
        rhs_cache.J_gas,
        rhs_cache.q_total,
        prop_cache.k,
        state,
        mesh,
        material,
    )

    # ========================================================================
    # Step 7: Compute divergences
    # ========================================================================
    divergence_1d!(rhs_cache.div_q, rhs_cache.q_total, state, mesh)
    divergence_species_1d!(rhs_cache.div_J, rhs_cache.J_gas, state, mesh)

    # ========================================================================
    # Step 8: Assemble final RHS
    # ========================================================================
    @inbounds for i = 1:n_cells_val
        # Energy equation: dT/dt = (-∇·q + S_conv + Q_rad + Q_rxn) / ρcp
        # (ρcp sums every component, pore gas included; the gas enthalpy flux
        # is carried by S_conv.)
        ρcp = prop_cache.ρcp[i]
        dT[i] =
            (
                -rhs_cache.div_q[i] +
                rhs_cache.S_conv[i] +
                rhs_cache.Q_rad[i] +
                rhs_cache.Q_rxn[i]
            ) / ρcp

        # Species equations - use unrolled version for type stability
        apply_species_residual_unrolled!(dξ, i, rhs_cache.div_J, rhs_cache.dξ_rxn, material)
    end

    return nothing
end

# ============================================================================
# Helper functions for boundary condition application
# ============================================================================

using ..BoundaryConditions: BoundaryConditionSet, get_thermal_bc, get_mass_bc
using ..Materials:
    RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY

"""
	apply_thermal_bc_flux!(q::Vector{Float64}, mesh,
						   bc_set::BoundaryConditionSet,
						   prop_cache,
						   material::Material{NC},
						   J_gas::Matrix,
						   t,
						   top_face_ids::Vector{Int},
						   bottom_face_ids::Vector{Int},
						   radiation_model::RadiationModel = SURFACE_ABSORPTION) where NC

Apply thermal boundary conditions using static dispatch (type-stable).

The surface energy balance solved here is `q_BC(T_s) = q_cond(T_s, T_cell)`.
Every `ConvectiveBC` in it applies its blowing model (`bc.blowing`, see
`blowing_state`) to the outward gas flux at the face, read from
`J_gas[g, face]` (which includes the Darcy venting part in Darcy–Fick mode).
`J_gas` is populated by the mass-flux pass that runs *before* this thermal BC
pass in `compute_rhs!`.

# Arguments
- `q`: Face flux array (modified in-place)
- `mesh`: The mesh
- `bc_set`: Type-stable boundary condition set
- `prop_cache`: Property cache
- `material`: Material (needed for per-component cp_g and is_gas predicate)
- `J_gas`: Gas mass-flux matrix `[NC × n_faces]`, with boundary entries populated
- `t`: Current time
- `top_face_ids`, `bottom_face_ids`: Pre-cached boundary face IDs
- `radiation_model`: Determines how `RadiativeFluxBC` handles absorption
"""
function apply_thermal_bc_flux!(
    q::AbstractVector,
    state,
    mesh,
    bc_set::BoundaryConditionSet,
    prop_cache,
    material::Material{NC},
    J_gas::AbstractMatrix,
    t,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int},
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
) where {NC}
    # Top boundary (boundary_id = 1, z = L)
    _apply_thermal_bc_static!(
        q,
        state,
        mesh,
        bc_set.thermal.top,
        prop_cache,
        material,
        J_gas,
        t,
        top_face_ids,
        Val(:top),
        radiation_model,
    )

    # Bottom boundary (boundary_id = 2, z = 0)
    _apply_thermal_bc_static!(
        q,
        state,
        mesh,
        bc_set.thermal.bottom,
        prop_cache,
        material,
        J_gas,
        t,
        bottom_face_ids,
        Val(:bottom),
        radiation_model,
    )

    return nothing
end


"""
Inner function for type-stable thermal BC application.

The `Val{boundary}` parameter enables compile-time dispatch on boundary type.
The `BC` type parameter enables compile-time dispatch on BC type.
The `radiation_model` determines how `RadiativeFluxBC` handles absorption.

For non-Dirichlet BCs, the boundary-temperature solver finds `T_s` from the
surface energy balance:
	q_BC(T_s) = k(T_s − T_cell)/Δz_half
"""
@inline function _apply_thermal_bc_static!(
    q::AbstractVector,
    state,
    mesh,
    bc::BC,
    prop_cache,
    material::Material{NC},
    J_gas::AbstractMatrix,
    t,
    face_ids::Vector{Int},
    ::Val{boundary},
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
) where {BC<:AbstractThermalBC,NC,boundary}
    # Handle DirichletThermalBC specially (T_s prescribed, no Newton solve)
    if BC <: DirichletThermalBC
        boundary_id = boundary === :top ? 1 : 2
        @inbounds for f in face_ids
            face = mesh.faces[f]
            left_id, right_id = face.cell_ids
            cell_id = left_id == 0 ? right_id : left_id
            T_cell = state.T[cell_id]

            T_bc = evaluate_param(bc.T, t)
            k_cell = prop_cache.k[cell_id]

            # Distance from cell center to boundary face (half cell width).
            # Uses the boundary face's own area — identical to the Newton
            # path's Δz_half = V/(A_face·2) — so the two paths agree under
            # a variable cross-section.
            dz_half = state.volume[cell_id] / state.A_face[f] / 2

            # Conductive flux: q = -k * dT/dz
            q_bc = k_cell * (T_bc - T_cell) / dz_half

            if boundary === :bottom
                q[f] = q_bc
            else
                q[f] = -q_bc
            end

            # Keep the diagnostic boundary state current (the bypass
            # previously skipped this, so surface-T diagnostics held
            # 300 K constructor placeholders for Dirichlet runs). Float64
            # storage, mirroring solve_boundary_temperature!.
            bstate = mesh.boundary_states[boundary_id]
            bstate.T_surface = ForwardDiff.value(T_bc)
            bstate.q_applied = ForwardDiff.value(q_bc)
            bstate.q_conducted = ForwardDiff.value(q_bc)
            bstate.blowing_mcp = 0.0
            bstate.blowing_ratio = 1.0
            bstate.pooling_factor = 1.0
            bstate.face_fraction = 0.0
            bstate.face_buoyancy = 1.0
            bstate.T_cell = ForwardDiff.value(T_cell)
            bstate.k_eff = ForwardDiff.value(k_cell)
            bstate.Δz_half = ForwardDiff.value(dz_half)
            bstate.converged = true
            bstate.iterations = 0
        end
    else
        # Boundary temperature solve for energy-conserving flux
        boundary_id = boundary === :top ? 1 : 2

        # Outward gas mass flux at the boundary face. `J_gas[:, f]` follows the
        # finite-volume +z convention, so escaping gas is +J at the top face and
        # -J at the bottom face.
        face_id_for_J = first(face_ids)
        J_outward = ntuple(NC) do j
            Jf = J_gas[j, face_id_for_J]
            boundary === :top ? Jf : -Jf
        end

        # Newton solve for T_s. The returned q_conducted is Tu-typed (carries
        # AD sensitivities through T_cell, k, J_gas, …).
        result = update_boundary_state!(
            state,
            mesh,
            boundary_id,
            bc,
            prop_cache,
            material,
            J_outward,
            t,
            radiation_model,
        )

        @inbounds for f in face_ids
            if boundary === :bottom
                q[f] = result.q_conducted
            else
                q[f] = -result.q_conducted
            end
        end
    end
    return nothing
end

"""
	apply_mass_bc_flux!(J::Matrix{Float64}, state, mesh, prop_cache,
						bc_set::BoundaryConditionSet, material::Material{NC}, t,
						top_face_ids::Vector{Int},
						bottom_face_ids::Vector{Int}) where NC

Apply mass boundary conditions using static dispatch (type-stable).

Only **gas** components receive a boundary flux: mass BCs model gas exchange
with the ambient (film law / venting), and a condensed phase cannot cross the
boundary. Writing the film law for solids (as an earlier revision did) put a
phantom O(100) kg/(m²·s) "solid boundary flux" into `J` that every production
consumer had to mask out — filtering at write time keeps `J`
physically meaningful for every reader, including diagnostics.

# Arguments
- `J`: Species flux matrix [NC × n_faces] (modified in-place)
- `state`: Geometry/state cache
- `mesh`: The mesh
- `prop_cache`: Property cache (for λ, ρ)
- `bc_set`: Type-stable boundary condition set
- `material`: Material (compile-time `is_gas` filter)
- `t`: Current time
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
"""
function apply_mass_bc_flux!(
    J::AbstractMatrix,
    state,
    mesh,
    prop_cache,
    bc_set::BoundaryConditionSet,
    material::Material{NC},
    t,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int},
) where {NC}
    # Top boundary (boundary_id = 1, z = L)
    _apply_mass_bc_static!(
        J,
        state,
        mesh,
        bc_set.mass.top,
        prop_cache,
        material,
        t,
        top_face_ids,
        Val(:top),
    )

    # Bottom boundary (boundary_id = 2, z = 0)
    _apply_mass_bc_static!(
        J,
        state,
        mesh,
        bc_set.mass.bottom,
        prop_cache,
        material,
        t,
        bottom_face_ids,
        Val(:bottom),
    )

    return nothing
end

"""
Inner function for type-stable mass BC application.
"""
@inline function _apply_mass_bc_static!(
    J::AbstractMatrix,
    state,
    mesh,
    bc::BC,
    prop_cache,
    material::Material{NC},
    t,
    face_ids::Vector{Int},
    ::Val{boundary},
) where {BC<:AbstractMassBC,boundary,NC}
    @inbounds for f in face_ids
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids
        cell_id = left_id == 0 ? right_id : left_id

        λ = prop_cache.λ[cell_id]
        ρ = prop_cache.ρ[cell_id]
        λρ = λ * ρ
        face_pos = state.face_positions[f]
        cell_ctr = state.cell_centers[cell_id]
        d = abs(face_pos - cell_ctr)

        # Sign convention: boundary-tag based. ACCUMULATE — in Darcy-Fick
        # mode the boundary face already carries the advective venting
        # flux from the pressure BC (compute_gas_fluxes_darcy_fick!);
        # the film-law diffusive part adds to it. In Fickian mode the
        # boundary faces arrive zero, so += reduces to a plain assignment.
        σ = boundary === :bottom ? 1.0 : -1.0
        _accumulate_mass_bc_gas_unrolled!(
            J,
            f,
            cell_id,
            bc,
            λρ,
            ρ,
            d,
            t,
            σ,
            state,
            material,
        )
    end
    return nothing
end

# Compile-time unrolled per-component accumulation; the `is_gas` predicate is
# resolved at codegen so solid/liquid rows are never touched.
@generated function _accumulate_mass_bc_gas_unrolled!(
    J::AbstractMatrix,
    f::Int,
    cell_id::Int,
    bc::AbstractMassBC,
    λρ::Real,
    ρ::Real,
    d::Real,
    t::Real,
    σ::Float64,
    state,
    material::Material{NC},
) where {NC}
    exprs = Expr[]
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                ξ_j = state.ξ[$j, cell_id]
                J[$j, f] += σ * mass_bc_flux(bc, ξ_j, λρ, ρ, d, t, $j)
            end
        end)
    end
    push!(exprs, :(return nothing))
    return Expr(:block, exprs...)
end


# ============================================================================
# Darcy-Fick RHS Computation
# ============================================================================

"""
	update_darcy_properties!(cache::PropertyCache, mesh::Mesh1D{NC},
							  material::Material{NC},
							  bc_set::BoundaryConditionSet,
							  top_face_ids::Vector{Int},
							  bottom_face_ids::Vector{Int}) where NC

Update Darcy-related properties using component-based permeability.

This method computes effective permeability from component κ values via the
material's `mixing.κ` rule. Components carry their own κ value (e.g.
`SolidComponent(:char, ..., κ=1e-12)`); a fully-solid material with κ=0 will
produce zero effective permeability and zero Darcy velocity.

The effective permeability is computed as a volume-weighted mixture of
component permeabilities, similar to thermal conductivity.

# Arguments
- `cache`: PropertyCache (modified in place)
- `mesh`: Mesh with current cell states
- `material`: Material definition with component κ values and `mixing.κ` spec
- `bc_set`: Type-stable BoundaryConditionSet containing mass BCs
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
"""
function update_darcy_properties!(
    cache::PropertyCache,
    state,
    mesh,
    material::Material{NC},
    bc_set::BoundaryConditionSet,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int},
) where {NC}
    n_cells = mesh.n_active_cells
    n_faces = n_cells + 1

    darcy = cache.darcy
    darcy === nothing && error(
        "update_darcy_properties!: PropertyCache has no Darcy sub-cache; build the workspace with a Darcy-enabled material.",
    )

    # Step 1: Compute cell-centered pressure, porosity, permeability, viscosity
    @inbounds for i = 1:n_cells
        T_i = state.T[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))

        # Porosity (gas volume fraction) — computed first so the ideal-gas pressure
        # can divide by it (gas occupies the pore space, not the bulk). Single
        # canonical implementation lives in `Properties.porosity`, which uses
        # intrinsic (skeletal) densities so the formula gives realistic values
        # for materials parameterized with bulk density + measured porosity.
        φ_i = porosity(ξ_i, T_i, material)
        darcy.φ[i] = φ_i

        # Pressure from ideal gas law in porous media: P = Σⱼ (ξⱼ R T / Mⱼ) / φ
        darcy.P[i] = Physics.ideal_gas_pressure(ξ_i, T_i, φ_i, material)

        # Permeability from component mixture (component-based)
        darcy.κ[i] = effective_permeability(ξ_i, T_i, material)

        # Pore-gas viscosity from the material's own curve (default: Sutherland air)
        darcy.μ[i] = Properties.mixture_gas_viscosity(ξ_i, T_i, material)
    end

    # Step 2: Compute face-centered pressure gradient and Darcy velocity
    @inbounds for f = 1:n_faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            # Boundary face: zero gradient (will be overwritten by BCs)
            darcy.∇P[f] = 0.0
            darcy.v_gas[f] = 0.0
            continue
        end

        # Distances from cell centers to face (distance-weighted harmonic mean)
        face_pos = state.face_positions[f]
        d_L = abs(face_pos - state.cell_centers[left_id])
        d_R = abs(state.cell_centers[right_id] - face_pos)
        d_LR = d_L + d_R

        # Pressure gradient
        darcy.∇P[f] = (darcy.P[right_id] - darcy.P[left_id]) / d_LR

        # Harmonic mean of mobility κ/μ at face
        κ_L = darcy.κ[left_id]
        κ_R = darcy.κ[right_id]
        μ_L = darcy.μ[left_id]
        μ_R = darcy.μ[right_id]

        # Mobility = κ/μ
        mob_L = μ_L > 1e-20 ? κ_L / μ_L : 0.0
        mob_R = μ_R > 1e-20 ? κ_R / μ_R : 0.0

        # Distance-weighted harmonic mean for series resistance (non-uniform meshes)
        # For uniform meshes (d_L = d_R), this reduces to standard harmonic mean.
        # An impermeable side is an infinite series resistance: gas cannot cross
        # into a κ ≈ 0 cell, so the face mobility is exactly 0.
        if mob_L > 1e-30 && mob_R > 1e-30
            mob_face = (d_L + d_R) / (d_L / mob_L + d_R / mob_R)
        else
            mob_face = 0.0
        end

        # Darcy velocity: v_g = -κ/μ × ∇P
        darcy.v_gas[f] = -mob_face * darcy.∇P[f]
    end

    # Step 3: Apply pressure boundary conditions using type-stable bc_set
    apply_pressure_bc_to_darcy_velocity!(
        darcy.v_gas,
        state,
        mesh,
        darcy.P,
        darcy.κ,
        darcy.μ,
        bc_set,
        top_face_ids,
        bottom_face_ids,
    )

    return cache
end

"""
	compute_rhs_darcy_fick!(dT::Vector{Float64}, dξ::Matrix{Float64},
							mesh::Mesh1D{NC}, prop_cache::PropertyCache,
							rhs_cache::RHSCache, material::Material{NC},
							bc_set::BoundaryConditionSet, t,
							top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int};
							I_radiation::Float64=0.0,
							radiation_model::RadiationModel=NO_RADIATION) where NC

Compute RHS with full Darcy-Fick gas transport using component-based permeability.

This version uses component κ values and mixing rules to compute effective
permeability, rather than a separate permeability model. Use this when:
- Components have κ values specified (e.g., SolidComponent(:char, ..., κ=1e-12))
- No legacy permeability_model is provided

# Arguments
- `dT`: Output vector of dT/dt [K/s]
- `dξ`: Output matrix of dξⱼ/dt [kg/(m³·s)], size [NC × n_cells]
- `mesh`: Mesh with current state
- `prop_cache`: Property cache (will be updated with Darcy properties)
- `rhs_cache`: RHS workspace cache
- `material`: Material definition with component κ values
- `bc_set`: Type-stable BoundaryConditionSet
- `t`: Current time [s]
- `top_face_ids`: Pre-cached face IDs at top boundary (boundary_id=1)
- `bottom_face_ids`: Pre-cached face IDs at bottom boundary (boundary_id=2)
- `I_radiation`: Surface radiation intensity [W/m²]
- `radiation_model`: Radiation model type

# Notes
This is the preferred API for Darcy-Fick transport. The effective permeability
is computed from component values via mixing rules (PARALLEL, SERIES, WEIGHTED).
Gas convective source is stored in rhs_cache.S_conv.
"""
function compute_rhs_darcy_fick!(
    dT::AbstractVector,
    dξ::AbstractMatrix,
    state,
    mesh::Mesh1D{NC},
    prop_cache::PropertyCache,
    rhs_cache::RHSCache,
    material::Material{NC},
    bc_set::BoundaryConditionSet,
    t,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int};
    I_radiation::Real = 0.0,
    radiation_model::RadiationModel = NO_RADIATION,
    radiation_kirchhoff::Bool = true,
) where {NC}

    n_cells_val = mesh.n_active_cells

    # ========================================================================
    # Step 1: Update Darcy properties using component-based permeability
    # ========================================================================
    update_darcy_properties!(
        prop_cache,
        state,
        mesh,
        material,
        bc_set,
        top_face_ids,
        bottom_face_ids,
    )

    # ========================================================================
    # Step 2: Compute conductive heat fluxes at interior faces
    # ========================================================================
    compute_conductive_fluxes!(rhs_cache.q_cond, state, mesh, prop_cache)
    copyto!(rhs_cache.q_total, rhs_cache.q_cond)

    # ========================================================================
    # Step 3: Reactions
    # ========================================================================
    @inbounds for i = 1:n_cells_val
        T = state.T[i]
        ξ = ntuple(j -> state.ξ[j, i], Val(NC))

        evaluate_reactions!(rhs_cache.rates, rhs_cache.heats, ξ, T, material)
        rhs_cache.Q_rxn[i] = reaction_heat_source(rhs_cache.rates, rhs_cache.heats)
        species_source_terms!(rhs_cache.species_src, rhs_cache.rates, material)
        for j = 1:NC
            rhs_cache.dξ_rxn[j, i] = rhs_cache.species_src[j]
        end
    end

    # ========================================================================
    # Step 4: Gas transport (Darcy + Fick) — interior + boundary
    # Must run BEFORE the thermal BC solve: the blowing models of the
    # convective BCs read the outward gas flux. Boundary faces receive the
    # advective venting flux from the pressure-BC Darcy velocity here; the
    # diffusive film-law part is then ADDED by `apply_mass_bc_flux!`.
    # ========================================================================
    Physics.compute_gas_fluxes_darcy_fick!(
        rhs_cache.J_gas,
        state,
        mesh,
        prop_cache,
        material,
        prop_cache.darcy.v_gas,
        prop_cache.darcy.φ,
    )

    apply_mass_bc_flux!(
        rhs_cache.J_gas,
        state,
        mesh,
        prop_cache,
        bc_set,
        material,
        t,
        top_face_ids,
        bottom_face_ids,
    )

    # ========================================================================
    # Step 5: Apply thermal boundary conditions
    # ========================================================================
    apply_thermal_bc_flux!(
        rhs_cache.q_total,
        state,
        mesh,
        bc_set,
        prop_cache,
        material,
        rhs_cache.J_gas,
        t,
        top_face_ids,
        bottom_face_ids,
        radiation_model,
    )

    # ========================================================================
    # Step 6: Compute radiation source (model-dependent)
    # ========================================================================
    compute_radiation_source!(
        rhs_cache.Q_rad,
        state,
        mesh,
        prop_cache,
        I_radiation,
        radiation_model,
        radiation_kirchhoff,
    )

    # ========================================================================
    # Step 7: Cell-centered gas convective source + divergences
    # (q_total carries the Step-5 thermal-BC boundary flux the boundary-cell
    # gradient needs.)
    # ========================================================================
    Physics.compute_gas_convective_source!(
        rhs_cache.S_conv,
        rhs_cache.J_gas,
        rhs_cache.q_total,
        prop_cache.k,
        state,
        mesh,
        material,
    )

    divergence_1d!(rhs_cache.div_q, rhs_cache.q_total, state, mesh)
    divergence_species_1d!(rhs_cache.div_J, rhs_cache.J_gas, state, mesh)

    # ========================================================================
    # Step 8: Assemble final RHS
    # ========================================================================
    @inbounds for i = 1:n_cells_val
        ρcp = prop_cache.ρcp[i]
        dT[i] =
            (
                -rhs_cache.div_q[i] +
                rhs_cache.S_conv[i] +
                rhs_cache.Q_rad[i] +
                rhs_cache.Q_rxn[i]
            ) / ρcp

        # Species equations
        apply_species_residual_unrolled!(dξ, i, rhs_cache.div_J, rhs_cache.dξ_rxn, material)
    end

    return nothing
end
