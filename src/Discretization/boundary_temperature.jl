# ============================================================================
# Boundary State Update Function
# ============================================================================

"""
	update_boundary_state!(state_cache, mesh, boundary_id, bc, cache, material,
	                       J_outward, t, radiation_model)

Solve for `T_s` satisfying the surface energy balance
`q_BC(T_s) = k_eff·(T_s − T_cell)/Δz_half` at boundary `boundary_id`
(1 = top, 2 = bottom). `J_outward[j]` is the per-component outward gas mass
flux at the boundary face (positive = leaving the slab), read by the blowing
models of the convective BCs. Auxiliary callers without the live `J_gas` pass
zeros, the no-blowing limit.
"""
function update_boundary_state!(
    state_cache,
    mesh::Mesh1D,
    boundary_id::Int,
    bc::AbstractThermalBC,
    cache::PropertyCache,
    material::Material{NC},
    J_outward::NTuple{NC,<:Real},
    t::Real,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
) where {NC}
    # Cell adjacent to this boundary (1 = top/exposed at z = L; 2 = bottom).
    if boundary_id == 1
        cell_id = mesh.n_active_cells
        face_id = mesh.cells[cell_id].face_ids[2]
    else
        cell_id = 1
        face_id = mesh.cells[cell_id].face_ids[1]
    end

    T_cell = state_cache.T[cell_id]
    k_eff = cache.k[cell_id]
    ε_eff = cache.ε[cell_id]
    a_eff = cache.a[cell_id]
    Δz_half = state_cache.volume[cell_id] / state_cache.A_face[face_id] / 2

    state = mesh.boundary_states[boundary_id]
    # `state` storage is Float64 (diagnostics); the returned tuple carries
    # Tu typing so AD derivatives flow through `q_conducted` into `q[f]`.
    return solve_boundary_temperature!(
        state,
        T_cell,
        k_eff,
        Δz_half,
        bc,
        material,
        J_outward,
        t,
        ε_eff;
        a_eff = a_eff,
        radiation_model = radiation_model,
    )
end

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
	check_boundary_energy_balance(mesh)

Verify that boundary energy balance is satisfied.

# Returns
Dict{Int, Float64}: Boundary ID to relative imbalance (should be ~machine precision).
"""
function check_boundary_energy_balance(mesh::AbstractMesh)
    results = Dict{Int,Float64}()

    for (boundary_id, state) in mesh.boundary_states
        if abs(state.q_applied) > 1e-10
            imbalance = abs(state.q_applied - state.q_conducted) / abs(state.q_applied)
        else
            imbalance = abs(state.q_applied - state.q_conducted)
        end
        results[boundary_id] = imbalance
    end

    return results
end

"""
	get_surface_temperature(mesh)

Get the top surface temperature from the solved boundary state.
"""
function get_surface_temperature(mesh::AbstractMesh)
    return mesh.boundary_states[1].T_surface
end

"""
	get_back_face_temperature(mesh)

Get the back face (bottom) temperature from the solved boundary state.
"""
function get_back_face_temperature(mesh::AbstractMesh)
    return mesh.boundary_states[2].T_surface
end
