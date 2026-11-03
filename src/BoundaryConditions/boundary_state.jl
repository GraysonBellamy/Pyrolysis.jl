# ============================================================================
# Boundary state (diagnostics of the surface energy balance)
# ============================================================================

"""
	BoundaryState

Stores solved boundary temperature and flux balance information.
Used to properly resolve the surface temperature that satisfies energy balance:
	q_applied(T_s) = q_conducted(T_s, T_cell)

# Fields
- `T_surface::Float64`: Solved surface temperature [K]
- `q_applied::Float64`: Heat flux from boundary condition [W/m²]
- `q_conducted::Float64`: Heat flux conducted into material [W/m²]
- `T_cell::Float64`: Adjacent cell temperature [K]
- `k_eff::Float64`: Effective conductivity at boundary [W/(m·K)]
- `Δz_half::Float64`: Distance from boundary to cell center [m]
- `converged::Bool`: Did Newton iteration converge?
- `iterations::Int`: Number of iterations used
- `blowing_mcp::Float64`: Blowing conductance `ṁ″·c_p,g = Σ_g ṁ_g·cp_g(T_s)`
  [W/(m²·K)] of the gas leaving the face at the converged `T_s`.
- `blowing_ratio::Float64`: Ratio `χ = h/h₀` the first convective member of the
  thermal BC applied at the converged `T_s` (1 without outflow or without a
  blowing model).
- `pooling_factor::Float64`, `face_fraction::Float64`, `face_buoyancy::Float64`:
  the pooling factor `ψ`, the face mass fraction of the volatiles and the face
  buoyancy `Λ` of a `PoolingFilm` model (1, 0 and 1 for every other model).
  See `BoundaryConditions.blowing_state`.
"""
mutable struct BoundaryState
    T_surface::Float64
    q_applied::Float64
    q_conducted::Float64
    T_cell::Float64
    k_eff::Float64
    Δz_half::Float64
    converged::Bool
    iterations::Int
    # Blowing diagnostics at the converged surface temperature: the blowing
    # conductance ṁ″·c_p,g [W/(m²·K)] of the gas leaving the face, the ratio
    # χ = h/h₀ the convective boundary applied, and for a pooling model the
    # pooling factor ψ, the face mass fraction of the volatiles and the face
    # buoyancy Λ (1, 0 and 1 for the other models).
    blowing_mcp::Float64
    blowing_ratio::Float64
    pooling_factor::Float64
    face_fraction::Float64
    face_buoyancy::Float64
end

# The surface-solve inputs and outcome; the blowing diagnostics start at their
# no-outflow values.
function BoundaryState(
    T_surface::Real,
    q_applied::Real,
    q_conducted::Real,
    T_cell::Real,
    k_eff::Real,
    Δz_half::Real,
    converged::Bool,
    iterations::Integer,
)
    return BoundaryState(
        T_surface,
        q_applied,
        q_conducted,
        T_cell,
        k_eff,
        Δz_half,
        converged,
        iterations,
        0.0,
        1.0,
        1.0,
        0.0,
        1.0,
    )
end

# Default constructor
BoundaryState() = BoundaryState(300.0, 0.0, 0.0, 300.0, 0.2, 1e-5, false, 0)
