"""
Boundary-face gas mass flux of the solver core.
"""

# ============================================================================
# Boundary-face gas mass flux (film law + Darcy venting)
# ============================================================================

"""
	_bc_boundary_face_fluxes(T_cell, ξ_cell, mass_bc, material,
	                         d_to_face, t, ::Val{tag}, darcy_active)
	    -> NTuple{NC}

Per-component boundary-face gas mass flux per area, **positive into the
domain**: the diffusive film law plus (when `darcy_active`) the pressure-BC
advective venting that mirrors the boundary-face branch of
`compute_gas_fluxes_darcy_fick!`. Non-gas components are zero. Generic over
Dual scalar types. Shared by the Jacobian's boundary operators
(`BCMassOperator`, `BCConvectionOperator`, `BCThermalOperator`) and the
solution extraction, so all of them see the flux the residual applies.
"""
function _bc_boundary_face_fluxes(
    T_cell::Real,
    ξ_cell::NTuple{NC,Tu},
    mass_bc,
    material::Material{NC},
    d_to_face::Real,
    t::Real,
    ::Val{tag},
    darcy_active::Bool,
) where {NC,Tu,tag}

    λ_eff = effective_gas_transfer(ξ_cell, T_cell, material)
    ρ_eff = mixture_density(ξ_cell, T_cell, material)
    λρ = λ_eff * ρ_eff

    Tv = promote_type(typeof(T_cell), Tu)
    F = ntuple(Val(NC)) do j
        if is_gas(material.components[j])
            Tv(mass_bc_flux(mass_bc, ξ_cell[j], λρ, ρ_eff, d_to_face, t, j))
        else
            zero(Tv)
        end
    end

    if darcy_active
        # Advective venting: outward Darcy speed advecting the intrinsic
        # density ξ/φ out of the boundary cell (inflow advects only the
        # zero-concentration ambient, hence the one-sided clamp).
        φ = porosity(ξ_cell, T_cell, material)
        P_int = ideal_gas_pressure(ξ_cell, T_cell, φ, material)
        κ_int = effective_permeability(ξ_cell, T_cell, material)
        μ_int = mixture_gas_viscosity(ξ_cell, T_cell, material)
        v_z = boundary_darcy_velocity(mass_bc, P_int, κ_int, μ_int, d_to_face, Val(tag))
        v_out = tag === :top ? v_z : -v_z
        if v_out > 0
            coef = v_out / max(φ, 1e-12)
            # Return a fresh tuple instead of rebinding `F` — a self-referential
            # rebind would box the variable (type-unstable, and inference can't
            # prove the boxed read defined).
            return ntuple(Val(NC)) do j
                if is_gas(material.components[j])
                    F[j] - coef * ξ_cell[j]
                else
                    F[j]
                end
            end
        end
    end

    return F
end
