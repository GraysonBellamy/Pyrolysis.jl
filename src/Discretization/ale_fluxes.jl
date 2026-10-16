"""
ALE mesh-motion advection of the temperature field.

The ALE time-derivative identity (chain rule for a moving sample point) is

    ∂f/∂t|_χ = ∂f/∂t|_x + w·∇f

where ∂/∂t|_x is the Eulerian (fixed-point) rate and ∂/∂t|_χ is the ALE
(mesh-point) rate.

Every component — condensed and gaseous — rides the skeleton, and the mesh
velocity is the skeleton velocity, so the species equations carry no
mesh-motion advection at all (their only mesh-motion term is the dilation
−ξθ applied in the residual). The thermal mass ρcₚ rides the mesh with them,
so the temperature is advected only by the velocity of the mesh RELATIVE to
the material, `w_rel = w_applied − w_physical`, which is identically zero in
production and non-zero only under the `Workspace.w_override` verification
hook. That mesh-relative term is what this file discretizes.

The w·∇T term is discretized with TVD-limited upwind differencing. Adding
w·∇T to the Eulerian RHS gives the mesh-frame equation ∂T/∂τ − w·∂T/∂z =
RHS — an advection equation with apparent transport velocity c = −w (the
lab-frame field drifts at −w relative to the moving mesh). The donor
(upwind) side therefore follows the sign of −w: for w > 0 the donor stencil
is {i, i+1, i+2} (right), for w < 0 it is {i-2, i-1, i} (left). A slope
limiter on the 3-point donor-side stencil achieves second-order accuracy on
smooth fields while preventing oscillations near discontinuities. The
limiter ψ(r) is applied to the near/far ratio of consecutive donor-side
gradients. Falls back to first-order at boundaries where the donor-donor
cell is unavailable. Default limiter: VANLEER.
"""

# Limiter, VANLEER, NONE, apply_limiter are available from flux_limiters.jl
# (included before this file in Discretization module)

"""
	compute_ale_advection_term!(advection::Vector{Float64},
								mesh::Mesh1D{NC},
								w_node::Vector{Float64}) where NC

Compute the ALE advection term w·∇T using upwind differencing.

From the ALE identity: ∂T/∂t|_χ = ∂T/∂t|_x + w·∇T
The result w·∇T is ADDED to the Eulerian RHS to form the ALE rate.

# Arguments
- `advection`: Output vector of (w·∇T)ᵢ [K/s] for each cell
- `mesh`: Mesh with current cell states (temperatures)
- `w_node`: Mesh velocities at nodes [m/s]

# Sign Convention
- w < 0: mesh contracting (nodes moving toward z=0)
- w > 0: mesh expanding (nodes moving away from z=0)
- The result is ADDED to the Eulerian dT/dt to get the ALE dT/dt

# Boundary Treatment
- Uses zero-gradient extrapolation (Neumann BC) at all boundaries

# Algorithm
Adding w·∇T to the Eulerian RHS gives the mesh-frame equation
∂T/∂τ − w·∂T/∂z = RHS: an advection equation with apparent transport
velocity c = −w (the lab-frame field drifts at −w relative to the mesh).
The donor (upwind) cell follows the sign of c, not the sign of w.
For cell i with mesh velocity wᵢ (average of bounding nodes):
- If wᵢ > 0 (c < 0): forward difference (T_right - T_i) / d  [donor on the right]
- If wᵢ < 0 (c > 0): backward difference (T_i - T_left) / d  [donor on the left]

In the standard shrinking case (w < 0, surface receding) the surface cell
therefore takes its gradient from the interior, not from the zero-gradient
boundary ghost.
"""
function compute_ale_advection_term!(
    advection::AbstractVector,
    state,
    mesh::Mesh1D{NC},
    w_node::AbstractVector;
    limiter::Limiter = VANLEER,
) where {NC}

    n_cells = mesh.n_active_cells

    @inbounds for i = 1:n_cells
        cell = mesh.cells[i]
        n_left, n_right = cell.node_ids

        w_cell = (w_node[n_left] + w_node[n_right]) / 2

        if abs(w_cell) < 1e-20
            advection[i] = zero(eltype(advection))
            continue
        end

        T_i = state.T[i]

        f_left, f_right = cell.face_ids
        face_left = mesh.faces[f_left]
        face_right = mesh.faces[f_right]

        left_id = face_left.cell_ids[1]
        right_id = face_right.cell_ids[2]

        T_left = _get_boundary_temperature(left_id, T_i, state)
        T_right = _get_boundary_temperature(right_id, T_i, state)

        if w_cell > 0
            # Apparent transport velocity c = -w < 0: donor is the right neighbor.
            right_right_id = _upwind_upwind_id(mesh, right_id, false)
            T_right_right = _get_boundary_temperature(right_right_id, T_right, state)
            d_i_R = _upwind_gradient_distance(state, i, right_id)
            d_R_RR =
                right_right_id != 0 ?
                abs(state.cell_centers[right_id] - state.cell_centers[right_right_id]) :
                d_i_R
            dT_dz = _limited_gradient(
                T_i,
                T_right,
                T_right_right,
                d_i_R,
                d_R_RR,
                right_right_id != 0,
                limiter,
                false,
            )
        else
            # Apparent transport velocity c = -w > 0: donor is the left neighbor.
            left_left_id = _upwind_upwind_id(mesh, left_id, true)
            T_left_left = _get_boundary_temperature(left_left_id, T_left, state)
            d_i_L = _upwind_gradient_distance(state, i, left_id)
            d_L_LL =
                left_left_id != 0 ?
                abs(state.cell_centers[left_id] - state.cell_centers[left_left_id]) : d_i_L
            dT_dz = _limited_gradient(
                T_i,
                T_left,
                T_left_left,
                d_i_L,
                d_L_LL,
                left_left_id != 0,
                limiter,
                true,
            )
        end

        advection[i] = w_cell * dT_dz
    end

    return nothing
end


"""
Helper function to get temperature at boundary for ALE advection.
Uses zero-gradient extrapolation (cell temperature) for boundary cells.

The zero-gradient approach is physically appropriate for the advection operator.
"""
@inline function _get_boundary_temperature(neighbor_id::Int, T_cell, state)
    if neighbor_id != 0
        return state.T[neighbor_id]
    end
    # Boundary face - use zero-gradient (cell temperature)
    return T_cell
end

"""
Distance across the upwind face, used as the denominator of the first-order
upwind gradient `(φ_i − φ_upwind)/d`. This is the state-cache analogue of the
mesh-level cell-center-to-face distance evaluated at the upwind face:
center-to-center for an interior neighbor, and the full cell thickness
`Δz_i = V_i/A_section` (= twice the center-to-face distance = the
center-to-mirror-ghost-cell spacing) for a boundary.

On non-uniform meshes this differs from a naive cell width and avoids the
zeroth-order `d/Δz` bias in the gradient (which would converge to (d/Δz)·f′
instead of f′ as h→0); reproducing that geometric distance keeps the boundary
stencil on the same center-to-ghost footing as the interior center-to-center
stencil.

At a boundary (neighbor_id == 0) the closed edge uses a zero-gradient ghost,
implemented by setting the upwind value equal to the cell's own value (see
`_get_boundary_temperature`). The first-order gradient is therefore
identically zero there, so this distance is immaterial to every consumer: in
the residual it divides a zero numerator. The full-width value is
nonetheless the correct, center-to-ghost-consistent one should the boundary
ever supply a non-trivial (inflow/Dirichlet) ghost.
"""
@inline function _upwind_gradient_distance(state, i::Int, neighbor_id::Int)
    if neighbor_id == 0
        # Boundary: center-to-ghost-cell distance = full cell thickness Δz_i
        # (= V_i/A_section = 2 × center-to-face). The returned value is inert in
        # practice — the zero-gradient ghost forces the gradient numerator to zero
        # (see above).
        return state.volume[i] / state.A_section[]
    end
    return abs(state.cell_centers[i] - state.cell_centers[neighbor_id])
end

"""
Get the upwind-upwind neighbor cell ID for second-order TVD reconstruction.
Returns 0 if the upwind-upwind cell doesn't exist (boundary or upwind is boundary).
"""
@inline function _upwind_upwind_id(mesh::Mesh1D, upwind_id::Int, upwind_from_left::Bool)
    if upwind_id == 0
        return 0
    end
    cell = mesh.cells[upwind_id]
    if upwind_from_left
        f = cell.face_ids[1]
        return mesh.faces[f].cell_ids[1]
    else
        f = cell.face_ids[2]
        return mesh.faces[f].cell_ids[2]
    end
end

"""
Compute TVD-limited spatial gradient ∂T/∂z at a cell center using donor-side
(upwind) reconstruction. The gradient is always in the +z direction regardless
of which side is the donor. The caller selects the donor side from the sign of
the apparent transport velocity (c = -w for the mesh-motion correction: w > 0
→ donor on the right, w < 0 → donor on the left).

Donor on the left: grad = (T_i - T_left) / d, positive when T increases with z
Donor on the right: grad = (T_right - T_i) / d, positive when T increases with z

The TVD correction uses a 3-point donor-side stencil:
	grad = grad_1 + β·ψ(r) × (grad_1 - grad_0),    β = d₀/(d₀ + d₁)

where grad_1 is the first-order donor-side gradient (in +z direction),
grad_0 is the same gradient one cell further to the donor side, d₀/d₁ are the
center-to-center distances of the near/far slope segments, and
r = grad_1/grad_0 is the near/far slope ratio. With this ratio the
correction stays bounded by the one-sided slopes on monotone data
(ψ(r) ≤ min(2r, 2) ⇒ |correction| ≤ |grad_1| terms only); the inverse ratio
grad_0/grad_1 would scale the correction with the far slope and is unbounded
as grad_0/grad_1 → ∞. On smooth fields r → 1, ψ(1) = 1 recovers the
second-order one-sided extrapolation grad_1 + β(grad_1 - grad_0): grad_1
estimates φ′ at the near-segment midpoint (d₀/2 to the donor side) and
(grad_1 − grad_0) ≈ ((d₀+d₁)/2)·φ″, so the distance-correct blend is
β = d₀/(d₀+d₁) — the uniform-mesh ½ constant would leave a first-order
¼(d₁−d₀)·φ″ error on stretched/AMR meshes.

Arguments use directional convention:
- T_i, T_upwind, T_upup: values at current, donor, and donor-donor cells
- upwind_from_left: if true, gradient = (T_i - T_upwind)/d; if false, = (T_upwind - T_i)/d

Falls back to first-order if the donor-donor cell doesn't exist or the far
slope is flat (|grad_0| < 1e-20, where the ratio is undefined).
"""
@inline function _limited_gradient(
    T_i::Real,
    T_upwind::Real,
    T_upup::Real,
    d_to_upwind::Real,
    d_upwind_to_upup::Real,
    has_upup::Bool,
    limiter::Limiter,
    upwind_from_left::Bool,
)

    # First-order gradient in +z direction
    if upwind_from_left
        # Upwind is to the left: grad = (T_i - T_left) / d
        grad_1 = (T_i - T_upwind) / d_to_upwind
    else
        # Upwind is to the right: grad = (T_right - T_i) / d
        grad_1 = (T_upwind - T_i) / d_to_upwind
    end

    if !has_upup || limiter == NONE
        return grad_1
    end

    # Previous gradient (one cell further upwind, same +z direction)
    if upwind_from_left
        grad_0 = (T_upwind - T_upup) / d_upwind_to_upup
    else
        grad_0 = (T_upup - T_upwind) / d_upwind_to_upup
    end

    if abs(grad_0) < 1e-20
        # Flat far-side slope: the limited correction is undefined (r → ±∞);
        # fall back to the bounded first-order gradient.
        return grad_1
    end
    r = grad_1 / grad_0

    ψ = apply_limiter(r, limiter)

    # Distance-correct blend: grad_1 sits d₀/2 to the donor side
    # of the cell center and (grad_1 − grad_0) ≈ ((d₀+d₁)/2)·φ″, so
    # β = d₀/(d₀+d₁) extrapolates exactly to the cell center (½ on uniform
    # meshes).
    β = d_to_upwind / (d_to_upwind + d_upwind_to_upup)
    return grad_1 + β * ψ * (grad_1 - grad_0)
end
