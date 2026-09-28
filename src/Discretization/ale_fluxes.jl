"""
ALE advection computation using upwind differencing.

In the ALE (Arbitrary Lagrangian-Eulerian) formulation, mesh motion is handled via
the ALE time derivative identity (chain rule for a moving sample point):

    ∂f/∂t|_χ = ∂f/∂t|_x + w·∇f

where ∂/∂t|_x is the Eulerian (fixed-point) rate and ∂/∂t|_χ is the ALE (mesh-point) rate.

The ALE equations are (Eulerian RHS + mesh correction):
- Energy: ∂T/∂t|_χ = (Eulerian RHS) + w·∇T      (T is intensive, no dilation)
- Gas species: ∂ξ/∂t|_χ = (Eulerian RHS) + w·∇ξ  (gas is spatially stationary)
- Mesh-bound species (solid/liquid): ∂ξ/∂t|_χ = S - ξ·θ  (dilation, handled separately)

The w·∇f advection terms are discretized using TVD-limited upwind differencing.
Adding w·∇f to the Eulerian RHS gives the mesh-frame equation
∂f/∂τ − w·∂f/∂z = RHS — an advection equation with apparent transport velocity
c = −w (the lab-frame field drifts at −w relative to the moving mesh). The donor
(upwind) side therefore follows the sign of −w: for w > 0 the donor stencil is
{i, i+1, i+2} (right), for w < 0 it is {i-2, i-1, i} (left). A slope limiter on
the 3-point donor-side stencil achieves second-order accuracy on smooth fields
while preventing oscillations near discontinuities. The limiter ψ(r) is applied
to the near/far ratio of consecutive donor-side gradients. Falls back to
first-order at boundaries where the donor-donor cell is unavailable.
Default limiter: VANLEER.
"""

using ..Materials: GAS, SOLID
# Limiter, VANLEER, NONE, apply_limiter are available from flux_limiters.jl
# (included before this file in Discretization module)

# ============================================================================
# ALE Convective Flux Computation
# ============================================================================

"""
	ale_convective_flux(φ_L::Float64, φ_R::Float64, v::Float64, w::Float64) -> Float64

Compute ALE convective flux at a face using first-order upwinding.

F = (v - w) * φ_upwind

# Arguments
- `φ_L, φ_R`: Field values in left and right cells
- `v`: Material velocity at face [m/s]
- `w`: Mesh velocity at face [m/s]

# Returns
Convective flux [φ units * m/s], positive = flow in +z direction
"""
@inline function ale_convective_flux(φ_L::Float64, φ_R::Float64, v::Float64, w::Float64)
    v_rel = v - w  # Relative velocity
    φ_face = v_rel >= 0.0 ? φ_L : φ_R
    return v_rel * φ_face
end

# ============================================================================
# Gas Velocity Computation (Darcy flow)
# ============================================================================

"""
	estimate_gas_velocity_from_generation!(v_gas::Vector{Float64}, mesh::Mesh1D{NC}, 
										   mass_generation_rate::Vector{Float64},
										   ρ_gas::Float64) where NC

Estimate gas velocity from mass generation (simplified model).

Integrate mass generation from bottom to surface to get velocity.
"""
function estimate_gas_velocity_from_generation!(
    v_gas::Vector{Float64},
    mesh::Mesh1D{NC},
    mass_generation_rate::Vector{Float64},
    ρ_gas::Float64,
) where {NC}
    n_faces_val = length(mesh.faces)
    n_cells_val = length(mesh.cells)

    # Integrate from bottom (last cell) toward surface
    # Mass flux at each face = integral of generation rate below
    cumulative_mass_flux = 0.0

    # Assume cells are ordered from surface (1) to bottom (n)
    # Gas flows toward surface (negative z direction in our convention)
    @inbounds for f = n_faces_val:-1:1
        face = mesh.faces[f]
        right_id = face.cell_ids[2]

        if right_id != 0
            # Add generation from this cell
            V = mesh.cell_volumes[right_id]
            A = mesh.face_areas[f]
            cumulative_mass_flux += mass_generation_rate[right_id] * V / A
        end

        # Velocity from mass flux
        if ρ_gas > 0.0
            v_gas[f] = -cumulative_mass_flux / ρ_gas  # Negative = toward surface
        else
            v_gas[f] = 0.0
        end
    end

    return nothing
end

# ============================================================================
# ALE Advection (Upwind Formulation)
# ============================================================================
#
# ALE time derivative identity (chain rule):
#     ∂f/∂t|_χ = ∂f/∂t|_x + w·∇f
#
# For the non-conservative temperature equation:
#     ρcₚ ∂T/∂t|_x = -∇·q + Q     (Eulerian)
#     ∂T/∂t|_χ = (-∇·q + Q)/ρcₚ + w·∇T     (ALE = Eulerian + w·∇T)
#
# The w·∇T term is ADDED to the Eulerian RHS (not subtracted). In standard
# form the mesh-frame equation is ∂T/∂τ + c·∂T/∂z = RHS with c = −w, so the
# upwind (donor) side follows the sign of −w, not the sign of w.
# ============================================================================

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
`_get_boundary_temperature` and the species/Darcy self-fill `col_upwind = i`).
The first-order gradient is therefore identically zero there, so this distance
is immaterial to every consumer: in the residual it divides a zero numerator,
and in the structured-Jacobian scatter the +w/d and −w/d face-local AD partials
fold into the same global column and cancel exactly (`_scatter_advection_row!`).
The full-width value is nonetheless the correct, center-to-ghost-consistent one
should the boundary ever supply a non-trivial (inflow/Dirichlet) ghost.
"""
@inline function _upwind_gradient_distance(state, i::Int, neighbor_id::Int)
    if neighbor_id == 0
        # Boundary: center-to-ghost-cell distance = full cell thickness Δz_i
        # (= V_i/A_section = 2 × center-to-face), matching the mesh boundary
        # face-distance convention. The returned value is inert in practice — the
        # zero-gradient ghost forces the gradient numerator to zero (see above).
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

"""
	compute_ale_species_advection!(advection, state, mesh::Mesh1D{NC},
								   material::Material{NC}, w_node;
								   θ_A = 0.0,
								   top_gas_open = false,
								   bottom_gas_open = false) where NC

Conservative flux-form mesh-motion transport for GAS components.

Gas is NOT mesh-bound (bulk gas velocity v_gas ≈ 0 in the lab frame), so a
moving mesh exchanges gas with its cells. The per-cell balance on the moving
volume is `d(ξ_g V)/dt = [ξ_g w A]_faces + (transport + sources)·V`, giving
the correction ADDED to the Eulerian RHS:

```
dξ_g/dt |mesh-motion = (ξ_face,R·w_R − ξ_face,L·w_L)/Δz − ξ_g,i·(V̇/V)
V̇/V = (w_R − w_L)/Δz + θ_A
```

The face-flux difference telescopes over the column, so **total gas mass
changes only through the boundary faces** — the discrete operator is exactly
conservative for any mesh motion. (A pointwise `+w·∇ξ_g` discretization would
not be: it does not telescope, so with shrinking cells and no escape path it
silently deletes the held gas with the vanishing volume.)

# Face values (first-order donor)
`ξ_face` is the donor-cell value on the side the face moves INTO
(`w_face > 0` → upper cell, `w_face < 0` → lower cell) — the gas the face
sweeps. First-order donor values keep the telescoping exact and match the
Darcy face-flux upwinding (`ale_gas_flux_darcy_fick`).

# Boundary faces
- **Impermeable mass BC** (`top_gas_open`/`bottom_gas_open` = false): gas
  cannot cross the face, so it rides the moving boundary — the face
  contributes NO sweep flux and the held gas compresses/dilutes through the
  `ξ·V̇/V` term (total mass exactly conserved in a fully closed column).
- **Open mass BC** (venting/flux/convective): zero-gradient ghost,
  `ξ_outside = ξ_cell` — a regressing surface releases the uncovered gas at
  the cell concentration; an advancing surface ingests at the same value
  (the mass BC's own diffusive/convective exchange handles far-field
  equilibration). This keeps the boundary discretization consistent with
  the continuum `∇·(ξw) − ξ∇·w = w·∇ξ`, which the MMS suite verifies.

Solid / liquid components are mesh-bound and receive NO advection term here;
their volume-change re-concentration `−ξⱼθ` is applied by the residual's
dilation step (tech-ref §11.2.3).

# Arguments
- `advection`: Output matrix [NC × n_cells] of correction terms [kg/(m³·s)]
- `state`: StateCache (reads ξ, volume, A_section)
- `mesh`: Mesh with connectivity
- `material`: Material definition (to determine solid vs gas components)
- `w_node`: Applied mesh velocities at nodes [m/s]
- `θ_A`: Column-global lateral strain rate (WithChi; 0 under NoChi) [1/s]
- `top_gas_open` / `bottom_gas_open`: whether the boundary mass BC lets gas
  cross (anything but `ImpermeableBC`)
"""
function compute_ale_species_advection!(
    advection::AbstractMatrix,
    state,
    mesh::Mesh1D{NC},
    material::Material{NC},
    w_node::AbstractVector;
    θ_A::Real = 0.0,
    top_gas_open::Bool = false,
    bottom_gas_open::Bool = false,
) where {NC}

    n_cells = mesh.n_active_cells
    ξ = state.ξ

    @inbounds for i = 1:n_cells
        cell = mesh.cells[i]
        n_left, n_right = cell.node_ids
        w_L = w_node[n_left]
        w_R = w_node[n_right]

        f_left, f_right = cell.face_ids
        left_id = mesh.faces[f_left].cell_ids[1]
        right_id = mesh.faces[f_right].cell_ids[2]

        if abs(w_L) < 1e-20 && abs(w_R) < 1e-20 && abs(θ_A) < 1e-20
            for j = 1:NC
                advection[j, i] = zero(eltype(advection))
            end
            continue
        end

        Δz = state.volume[i] / state.A_section[]

        _compute_species_fluxform_unrolled!(
            advection,
            ξ,
            i,
            left_id,
            right_id,
            w_L,
            w_R,
            Δz,
            θ_A,
            top_gas_open,
            bottom_gas_open,
            material,
        )
    end

    return nothing
end


# Conservative flux-form gas transport (compile-time unrolled for type
# stability). Reads ξ values directly from the state matrix using explicit
# cell indices, avoiding closure allocations.
@generated function _compute_species_fluxform_unrolled!(
    advection::AbstractMatrix,
    ξ::AbstractMatrix,
    i::Int,
    left_id::Int,
    right_id::Int,
    w_L::Real,
    w_R::Real,
    Δz::Real,
    θ_A::Real,
    top_gas_open::Bool,
    bottom_gas_open::Bool,
    material::Material{NC},
) where {NC}
    exprs = Expr[]
    for j = 1:NC
        push!(
            exprs,
            quote
                @inbounds if is_gas(material.components[$j])
                    ξ_i = ξ[$j, i]

                    # Sweep flux ξ_face·w through each face; donor = side the
                    # face moves into. Boundary faces per the docstring.
                    flux_R = if right_id != 0
                        (w_R > 0 ? ξ[$j, right_id] : ξ_i) * w_R
                    elseif top_gas_open
                        ξ_i * w_R           # zero-gradient ghost
                    else
                        zero(ξ_i * w_R)     # impermeable: gas rides the face
                    end
                    flux_L = if left_id != 0
                        (w_L > 0 ? ξ_i : ξ[$j, left_id]) * w_L
                    elseif bottom_gas_open
                        ξ_i * w_L
                    else
                        zero(ξ_i * w_L)
                    end

                    advection[$j, i] =
                        (flux_R - flux_L) / Δz - ξ_i * ((w_R - w_L) / Δz + θ_A)
                else
                    # SOLID and LIQUID components move with the mesh (v = w):
                    # no relative flux; their dilation term lives in the
                    # residual's dedicated correction step.
                    advection[$j, i] = zero(eltype(advection))
                end
            end,
        )
    end
    return Expr(:block, exprs...)
end

# ============================================================================
# ALE Darcy-Fick Gas Fluxes
# ============================================================================
#
# In the ALE frame, gas advects at the RELATIVE INTERSTITIAL velocity. The
# lab-frame advective flux is (ξ/φ)·v_gas (ξ is per bulk volume, the Darcy
# velocity is superficial), and the moving face sweeps bulk concentration ξ·w,
# so the relative flux is exactly
#     (ξ/φ)·v_gas − ξ·w = ξ · (v_gas/φ − w)
#
# The combined Darcy-Fick flux in ALE frame is:
#     Jⱼ = ξⱼ × (v_gas/φ_face - w) - (λ_face/T_face) × ∂(ξⱼT)/∂x
#          ↑ advection (ALE)         ↑ diffusion (volume fraction gradient)
#
# This properly accounts for:
# - Gas moving through the porous medium (Darcy flow, v_gas/φ interstitial)
# - Mesh motion (ALE, w)
# - Diffusive transport (volume fraction gradient)
# ============================================================================

"""
	mesh_velocity_at_face(w_node::Vector{Float64}, mesh::Mesh1D{NC}, f::Int) -> Float64

Get mesh velocity at face f.

In 1D, the face is located at a single node, so we use the node velocity directly.
"""
@inline function mesh_velocity_at_face(
    w_node::Vector{Float64},
    mesh::Mesh1D{NC},
    f::Int,
) where {NC}
    node_id = mesh.faces[f].node_id
    return w_node[node_id]
end

"""
	ale_gas_flux_darcy_fick(v_g::Float64, w::Float64, ξ_L::Float64, ξ_R::Float64,
							φ_L::Float64, φ_R::Float64,
							λ_face::Float64, T_L::Float64, T_R::Float64,
							d_LR::Float64) -> Float64

Compute gas flux at a face in ALE frame using Darcy-Fick formulation [kg/(m²·s)].

Combines advective transport at the relative interstitial velocity with
diffusive transport:
	J = ξ_upwind × (v_g/φ_face - w) - (λ_face/T_face) × (ξ_R×T_R - ξ_L×T_L) / d_LR
		↑ advection (ALE)             ↑ diffusion (volume fraction gradient)

This is the exact identity (ξ/φ)·v_g − ξ·w = ξ·(v_g/φ − w): the superficial
Darcy velocity advects the intrinsic density ξ/φ while the moving face sweeps
bulk concentration. φ_face = (φ_L + φ_R)/2, floored at 1e-12 consistently with
the pressure closure.

# Arguments
- `v_g`: superficial Darcy gas velocity at face [m/s]
- `w`: Mesh velocity at face [m/s]
- `ξ_L, ξ_R`: Mass concentrations on left and right [kg/m³]
- `φ_L, φ_R`: Porosity on left and right [-]
- `λ_face`: Face transport coefficient (harmonic mean of λ) [m²/s]
- `T_L, T_R`: Temperatures on left and right [K]
- `d_LR`: Distance between cell centers [m]

# Returns
Total flux [kg/(m²·s)], positive = flow in +z direction
"""
@inline function ale_gas_flux_darcy_fick(
    v_g::Float64,
    w::Float64,
    ξ_L::Float64,
    ξ_R::Float64,
    φ_L::Float64,
    φ_R::Float64,
    λ_face::Float64,
    T_L::Float64,
    T_R::Float64,
    d_LR::Float64,
)
    # Relative interstitial velocity in ALE frame
    φ_face = max(0.5 * (φ_L + φ_R), 1e-12)
    v_rel = v_g / φ_face - w

    # Diffusive flux (volume fraction / partial pressure gradient)
    T_face = 0.5 * (T_L + T_R)
    J_diff = -λ_face * (ξ_R * T_R - ξ_L * T_L) / (T_face * d_LR)

    # Advective flux with first-order upwinding based on relative velocity
    ξ_upwind = v_rel >= 0.0 ? ξ_L : ξ_R
    J_adv = v_rel * ξ_upwind

    return J_adv + J_diff
end

"""
	compute_ale_gas_fluxes_darcy_fick!(J::Matrix{Float64}, mesh::Mesh1D{NC},
									   cache::PropertyCache, material::Material{NC},
									   v_gas::Vector{Float64}, φ::Vector{Float64},
									   w_node::Vector{Float64}) where NC

Compute gas mass fluxes in ALE frame using full Darcy-Fick formulation.

This combines:
- Darcy advective transport (gas flow through porous medium)
- ALE mesh motion correction
- Diffusive transport (volume fraction gradient)

The flux formula is:
	Jⱼ = ξⱼ × (v_gas/φ_face - w) - (λ_face/T_face) × ∂(ξⱼT)/∂x

with φ_face the arithmetic-mean face porosity — see
[`ale_gas_flux_darcy_fick`](@ref) for the derivation.

# Arguments
- `J`: Output matrix J[j,f] = flux of component j at face f [kg/(m²·s)]
- `mesh`: Mesh with current state
- `cache`: Property cache (must have λ_face computed)
- `material`: Material definition
- `v_gas`: superficial Darcy gas velocity at faces [m/s]
- `φ`: Porosity at cell centers [-] (e.g. `cache.darcy.φ`)
- `w_node`: Mesh velocity at nodes [m/s]

# Notes
- When w_node = 0, reduces to standard Darcy-Fick
- When v_gas = 0, gives ALE advection due to mesh motion only
- Solid components have zero flux (transport handled by mesh motion)
"""
function compute_ale_gas_fluxes_darcy_fick!(
    J::AbstractMatrix,
    state,
    mesh::Mesh1D{NC},
    cache::PropertyCache,
    material::Material{NC},
    v_gas::AbstractVector,
    φ::AbstractVector,
    w_node::AbstractVector,
) where {NC}
    n_faces_val = mesh.n_active_cells + 1
    ξ = state.ξ

    fill!(J, zero(eltype(J)))

    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            continue
        end

        w_face = mesh_velocity_at_face(w_node, mesh, f)

        d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])
        λ_face = cache.λ_face[f]
        T_L = state.T[left_id]
        T_R = state.T[right_id]
        v_g = v_gas[f]
        φ_L = φ[left_id]
        φ_R = φ[right_id]

        _compute_ale_gas_fluxes_darcy_fick_unrolled!(
            J,
            ξ,
            f,
            left_id,
            right_id,
            v_g,
            w_face,
            φ_L,
            φ_R,
            λ_face,
            T_L,
            T_R,
            d_LR,
            material,
        )
    end

    return nothing
end


"""
Compile-time unrolled ALE Darcy-Fick flux computation for type stability.
Reads ξ values directly from the state matrix to avoid closure allocations.
"""
@generated function _compute_ale_gas_fluxes_darcy_fick_unrolled!(
    J::Matrix,
    ξ::AbstractMatrix,
    f::Int,
    left_id::Int,
    right_id::Int,
    v_g,
    w_face,
    φ_L,
    φ_R,
    λ_face,
    T_L,
    T_R,
    d_LR,
    material::Material{NC},
) where {NC}

    exprs = Expr[quote
        # Relative interstitial velocity in ALE frame: the superficial Darcy
        # velocity advects intrinsic density ξ/φ; the moving face sweeps bulk
        # ξ — exactly ξ·(v_g/φ_face − w).
        φ_face = max(0.5 * (φ_L + φ_R), 1e-12)
        v_rel = v_g / φ_face - w_face
    end]
    for j = 1:NC
        push!(
            exprs,
            quote
                @inbounds if is_gas(material.components[$j])
                    ξ_L_j = ξ[$j, left_id]
                    ξ_R_j = ξ[$j, right_id]

                    # Diffusive flux (volume fraction gradient)
                    T_face = 0.5 * (T_L + T_R)
                    J_diff = -λ_face * (ξ_R_j * T_R - ξ_L_j * T_L) / (T_face * d_LR)

                    # Advective flux with relative velocity (first-order upwind)
                    ξ_upwind = v_rel >= 0.0 ? ξ_L_j : ξ_R_j
                    J_adv = v_rel * ξ_upwind

                    # Combined ALE Darcy-Fick flux
                    J[$j, f] = J_adv + J_diff
                end
            end,
        )
    end
    return Expr(:block, exprs...)
end
