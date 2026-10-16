"""
P1 approximation for volumetric radiation transport.

Implements the first-order spherical harmonics (P1) method for solving
radiative transfer in optically thick to intermediate media. This is the
industry-standard approach for pyrolysis modeling where optical thickness
varies during decomposition.

The P1 method solves for incident radiation intensity G [W/m²]:
	∇·(D_rad ∇G) = α G - 4σαT⁴

where:
- D_rad = 1/(3α) is the radiation diffusion coefficient
- α is absorption coefficient [m⁻¹]
- ε is emissivity [-]
- σ is Stefan-Boltzmann constant

Reference: 
- P1_Volumetric_Reradiation_Implementation.md
- Modest, M.F., Radiative Heat Transfer, 3rd Ed., 2013, Chapter 15
"""

# ============================================================================
# Core P1 Functions
# ============================================================================

# Transparent-limit floor on the extinction coefficient [m⁻¹]. With D = 0 the
# quasi-steady Newton system loses its G row entirely (singular tridiagonal
# matrix where α = 0), so the limit is taken as a large-but-finite diffusion
# coefficient instead. This is the single canonical guard — every consumer of
# D_rad goes through `compute_radiation_diffusion_coefficient`.
const P1_MIN_EXTINCTION = 1e-10

"""
	compute_radiation_diffusion_coefficient(α::Float64, σ_s::Float64=0.0) -> Float64

Compute the P1 radiation diffusion coefficient.

	D_rad = 1 / (3 · max(α + σ_s, $(P1_MIN_EXTINCTION)))

For pyrolysis materials, scattering is typically negligible (σ_s ≈ 0):
	D_rad = 1 / (3α)

# Arguments
- `α`: Absorption coefficient [m⁻¹].
- `σ_s`: Scattering coefficient [m⁻¹] (default: 0). Typically zero for pyrolysis.

# Returns
- `D_rad`: Radiation diffusion coefficient [m]

# Physical Interpretation
D_rad represents the "mean free path" for radiation diffusion. Large D_rad
(small α) means radiation travels far before being absorbed. Small D_rad
(large α) means radiation is absorbed locally.

# Transparent limit
Below an extinction of `P1_MIN_EXTINCTION` the coefficient is clamped to the
large finite value `1/(3·P1_MIN_EXTINCTION)` rather than diverging (or being
zeroed): a transparent cell then behaves as a near-perfect conductor of G,
which keeps the quasi-steady Newton system regular. P1 itself is not a valid
model in this regime (τ < 0.1) — see [`check_p1_validity`](@ref).

# Examples
```julia
# Typical virgin polymer (α ≈ 1000-2000 m⁻¹)
D_rad = compute_radiation_diffusion_coefficient(1500.0)  # ≈ 2.2×10⁻⁴ m

# Char layer (α ≈ 5000-10000 m⁻¹)
D_rad = compute_radiation_diffusion_coefficient(8000.0)  # ≈ 4.2×10⁻⁵ m

# With scattering (rare for pyrolysis)
D_rad = compute_radiation_diffusion_coefficient(5000.0, 1000.0)
```
"""
@inline function compute_radiation_diffusion_coefficient(α::Float64, σ_s::Float64 = 0.0)
    return 1.0 / (3.0 * max(α + σ_s, P1_MIN_EXTINCTION))
end

"""
	compute_radiation_source_terms!(Q_emission::Vector{Float64},
									Q_absorption::Vector{Float64},
									mesh::Mesh1D{NC},
									rad_cache::RadiationCache,
									prop_cache::PropertyCache) where NC

Compute volumetric emission and absorption source terms for P1 equation.

# Emission Source (blackbody emission into radiation field):
	Q_emission[i] = 4σ × ε[i] × α[i] × T[i]⁴    [W/m³]

# Absorption (removal from radiation field):
	Q_absorption[i] = α[i] × G[i]                [W/m³]

The net radiation source for the energy equation is:
	Q_net = Q_absorption - Q_emission

Note: These are the sources for the G equation. For the T equation,
the radiation contribution comes through ∇·(D_rad ∇G).

# Arguments
- `Q_emission`: Output vector for emission sources [W/m³]
- `Q_absorption`: Output vector for absorption sources [W/m³]
- `mesh`: Mesh with current temperature state
- `rad_cache`: Radiation cache with G values
- `prop_cache`: Property cache with α, ε values
"""
function compute_radiation_source_terms!(
    Q_emission::Vector{Float64},
    Q_absorption::Vector{Float64},
    mesh,
    rad_cache::RadiationCache,
    prop_cache::PropertyCache,
)
    n_cells = length(mesh.cells)

    @inbounds for i = 1:n_cells
        T = mesh.cell_states[i].T
        G = rad_cache.G[i]
        # Use volumetric α from rad_cache if populated, otherwise from prop_cache
        # Note: prop_cache.α is already volumetric (effective_absorption = Σⱼ ξⱼ αⱼ [m⁻¹])
        α =
            length(rad_cache.α) >= i && rad_cache.α[i] > 0 ? rad_cache.α[i] :
            prop_cache.α[i]

        # Emission: 4σαT⁴ (blackbody emission, NO emissivity factor)
        # Physical: Energy emitted per unit volume into radiation field
        # Note: In standard P1 approximation, volumetric emission is blackbody.
        # Emissivity only appears at boundaries (Marshak BC) where it represents
        # surface reflectivity/absorptivity. The interior is assumed to be a
        # perfect absorber/emitter (gray body with ε=1 locally).
        Q_emission[i] = 4.0 * STEFAN_BOLTZMANN * α * T^4

        # Absorption: αG
        # Physical: Energy absorbed from radiation field per unit volume
        Q_absorption[i] = α * G
    end

    return nothing
end

"""
	update_radiation_properties!(rad_cache::RadiationCache, mesh,
								 prop_cache::PropertyCache)

Update radiation-specific properties in the cache. This is the **single**
method of this function — the transparent limit is handled exclusively by
[`compute_radiation_diffusion_coefficient`](@ref) (large finite D, never 0).

Computes:
1. Cell-centered D_rad[i] = 1/(3α[i])
2. Face-centered D_rad_face[f] = harmonic mean weighted by distance

The harmonic mean is used for D_rad at faces because radiation transport
follows resistance-in-series (like conduction). This ensures correct
flux at interfaces between materials with different optical properties.

# Arguments
- `rad_cache`: Radiation cache to update
- `mesh`: Mesh with geometry information
- `prop_cache`: Property cache with α values (must be updated first)
"""
function update_radiation_properties!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
)
    n_cells = mesh.n_active_cells
    n_faces = n_cells + 1

    # Ensure cache is properly sized
    if length(rad_cache.D_rad) != n_cells
        resize_radiation_cache!(rad_cache, n_cells)
    end

    # Cell-centered D_rad and volumetric α
    @inbounds for i = 1:n_cells
        # prop_cache.α is already volumetric [m⁻¹] from effective_absorption = Σⱼ ξⱼ αⱼ
        # where ξⱼ [kg/m³] × αⱼ [m²/kg] = [m⁻¹]
        α_vol = prop_cache.α[i]

        # Store volumetric α in rad_cache
        rad_cache.α[i] = α_vol

        # D_rad = 1/(3α), transparent-limit clamped in the canonical helper
        rad_cache.D_rad[i] = compute_radiation_diffusion_coefficient(α_vol)
    end

    # Face-centered D_rad using distance-weighted harmonic mean
    @inbounds for f = 1:n_faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids
        face_pos = mesh.face_positions[f]

        if left_id == 0
            # Left boundary: use right cell value
            rad_cache.D_rad_face[f] = rad_cache.D_rad[right_id]
        elseif right_id == 0
            # Right boundary: use left cell value
            rad_cache.D_rad_face[f] = rad_cache.D_rad[left_id]
        else
            # Interior face: distance-weighted harmonic mean
            d_L = abs(face_pos - mesh.cell_centers[left_id])
            d_R = abs(mesh.cell_centers[right_id] - face_pos)

            rad_cache.D_rad_face[f] = harmonic_mean_weighted(
                rad_cache.D_rad[left_id],
                rad_cache.D_rad[right_id],
                d_L,
                d_R,
            )
        end
    end

    return nothing
end

# ============================================================================
# Radiation Flux Computation
# ============================================================================

"""
	compute_radiation_fluxes!(rad_cache::RadiationCache,
							  mesh::Mesh1D{NC},
							  prop_cache::PropertyCache,
							  radiation_bcs::RadiationBCSet,
							  t) where NC

Compute radiation diffusion fluxes at all faces using P1 approximation.

For interior faces:
	q_rad[f] = -D_rad_face × (G[R] - G[L]) / Δx

For boundary faces, the Marshak or other BC is applied.

# Arguments
- `rad_cache`: Radiation cache (updated in-place: q_rad)
- `mesh`: Mesh with geometry
- `prop_cache`: Property cache with α, ε values
- `radiation_bcs`: RadiationBCSet with radiation boundary conditions
- `t`: Current simulation time [s] (supports dual numbers for AD)

# Notes
- Type-generic to support automatic differentiation
"""
function compute_radiation_fluxes!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
    radiation_bcs::RadiationBCSet,
    t,
)
    n_faces = length(mesh.faces)

    # Interior faces (radiation diffusion)
    @inbounds for f = 1:n_faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            continue
        end

        D_rad_face = rad_cache.D_rad_face[f]
        Δx = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
        rad_cache.q_rad[f] =
            -D_rad_face * (rad_cache.G[right_id] - rad_cache.G[left_id]) / Δx
    end

    _apply_all_radiation_bcs!(rad_cache, mesh, prop_cache, radiation_bcs, t)

    return nothing
end

# RadiationBCSet (type-stable) boundary application
@inline function _apply_all_radiation_bcs!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
    radiation_bcs::RadiationBCSet,
    t,
)
    faces_top = get(mesh.boundary_faces, 1, nothing)
    if faces_top !== nothing && !isempty(faces_top)
        _apply_bc!(
            rad_cache,
            mesh,
            prop_cache,
            radiation_bcs.top,
            faces_top::Vector{Int},
            t,
        )
    end

    faces_bottom = get(mesh.boundary_faces, 2, nothing)
    if faces_bottom !== nothing && !isempty(faces_bottom)
        _apply_bc!(
            rad_cache,
            mesh,
            prop_cache,
            radiation_bcs.bottom,
            faces_bottom::Vector{Int},
            t,
        )
    end
    return nothing
end

# Dispatch helpers for RadiationBCSet (method-based, no isa chain)
@inline function _apply_bc!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
    bc::MarshakRadiationBC,
    faces::Vector{Int},
    t,
)
    _apply_marshak_bc!(rad_cache, mesh, prop_cache, bc, faces, t)
    return nothing
end

@inline function _apply_bc!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
    bc::ReflectiveRadiationBC,
    faces::Vector{Int},
    t,
)
    _apply_reflective_bc!(rad_cache, faces)
    return nothing
end

@inline function _apply_bc!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
    bc::DirichletRadiationBC,
    faces::Vector{Int},
    t,
)
    _apply_dirichlet_bc!(rad_cache, faces)
    return nothing
end

# ============================================================================
# Marshak BC kernels — canonical coefficient / driving-field helpers
# ============================================================================
#
# Derivation (half-range fluxes of the P1 intensity, inward normal n̂):
#   q_in  = G/4 + q/2,   q_out = G/4 − q/2          (q = net flux along n̂)
# A gray diffuse surface of emissivity ε transmits ε of the external radiosity
# (I_ext + σT_amb⁴) and reflects (1−ε) of the outgoing half-range flux:
#   q_in = ε(I_ext + σT_amb⁴) + (1−ε)·q_out
# Solving for the net flux into the domain:
#   q = c · (G_ext + 4·I_ext − G),   c = ε/(2(2−ε)),   G_ext = 4σT_amb⁴
# The external incident flux therefore carries the weight 4c = 2ε/(2−ε)
# (2.0 at ε = 1, 1.64 at ε = 0.9) — the same weight the ambient term carries
# (c·4σT_amb⁴ ≡ 4c·σT_amb⁴). It is NOT 1.

"""
	_marshak_exchange_coefficient(bc::MarshakRadiationBC, T_surf) -> c

The Marshak G-exchange coefficient c = ε/(2(2−ε)) (Modest, Radiative Heat
Transfer, 3rd ed., Eq. 16.48), with ε evaluated at `T_surf`.
"""
@inline function _marshak_exchange_coefficient(bc::MarshakRadiationBC, T_surf)
    ε = evaluate_param_T(bc.emissivity, T_surf)
    return ε / (2.0 * (2.0 - ε))
end

"""
	_marshak_driving_field(bc::MarshakRadiationBC, t) -> G_ext + 4·I_ext

The driving incident field of the Marshak relation
`q = c·(G_ext + 4·I_ext − G)`: the ambient equilibrium field
G_ext = 4σT_amb⁴ plus the external incident flux at its half-range weight 4
(so that c·4·I_ext = 2ε/(2−ε)·I_ext).
"""
@inline function _marshak_driving_field(bc::MarshakRadiationBC, t)
    I_ext = evaluate_param(bc.I_external, t)
    T_amb = evaluate_param(bc.T_ambient, t)
    return 4.0 * STEFAN_BOLTZMANN * T_amb^4 + 4.0 * I_ext
end

"""
	_marshak_effective_coefficient(bc, T_cell, D_cell, δ) -> c_eff

Effective Marshak exchange coefficient seen by the **cell-centered** value
G_N of the boundary cell: the surface exchange conductance c = ε/(2(2−ε)) in
series with the half-cell diffusive conductance β = D_rad/δ (δ = distance
from the cell center to the boundary face):

	c_eff = c·β/(c + β) = 1 / (1/c + δ/D_rad)

so the applied flux `q = c_eff·(G_ext + 4·I_ext − G_N)` is algebraically
identical to the exact face relation `q = c·(G_ext + 4·I_ext − G_face)`
closed with the discrete diffusive flux q = D_rad·(G_face − G_N)/δ. This
removes the O(δ) error of evaluating the Marshak relation at the cell
center. Used by both the flux application (`_apply_marshak_bc!`) and the
quasi-steady Newton Jacobian — the two must stay consistent.
"""
@inline function _marshak_effective_coefficient(bc::MarshakRadiationBC, T_cell, D_cell, δ)
    c = _marshak_exchange_coefficient(bc, T_cell)
    β = D_cell / max(δ, eps())
    return c * β / (c + β)
end

# Type-specific implementations - these are called with concrete types
@inline function _apply_marshak_bc!(
    rad_cache::RadiationCache,
    mesh,
    prop_cache::PropertyCache,
    bc::MarshakRadiationBC,
    faces::Vector{Int},
    t,
)
    @inbounds for f in faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids
        cell_id = left_id == 0 ? right_id : left_id

        T_cell = mesh.cell_states[cell_id].T
        G_cell = rad_cache.G[cell_id]
        D_cell = rad_cache.D_rad[cell_id]
        δ = abs(mesh.face_positions[f] - mesh.cell_centers[cell_id])

        # Half-cell-corrected Marshak flux in terms of the cell-centered G:
        # q = c_eff·(G_ext + 4·I_ext − G_cell). ε(T) is evaluated at the
        # boundary-cell temperature (the surface-Newton T_s is not available
        # inside the radiation solve).
        c_eff = _marshak_effective_coefficient(bc, T_cell, D_cell, δ)
        flux_into = c_eff * (_marshak_driving_field(bc, t) - G_cell)
        if right_id == 0
            rad_cache.q_rad[f] = -flux_into
        else
            rad_cache.q_rad[f] = flux_into
        end
    end
    return nothing
end

@inline function _apply_reflective_bc!(rad_cache::RadiationCache, faces::Vector{Int})
    # Reflective BC: zero flux
    @inbounds for f in faces
        rad_cache.q_rad[f] = 0.0
    end
    return nothing
end

@inline function _apply_dirichlet_bc!(rad_cache::RadiationCache, faces::Vector{Int})
    # Dirichlet BC: flux determined by interior gradient (handled in discretization)
    # For now, set to zero as placeholder
    @inbounds for f in faces
        rad_cache.q_rad[f] = 0.0
    end
    return nothing
end

# ============================================================================
# Boundary Condition Flux Functions
# ============================================================================

"""
	radiation_bc_flux(bc::MarshakRadiationBC, G_surf, T_surf, α_surf, ε_surf, t)

Compute radiation flux at boundary with Marshak condition, evaluated at the
**surface (face) value** of G.

The P1 Marshak condition for a gray diffuse surface with emissivity ε:

	q = c · (G_ext + 4·I_ext - G_surf),    c = ε / (2(2-ε))

where:
- I_ext is external incident radiation [W/m²], carrying its half-range
  weight 4c = 2ε/(2-ε) — the **same** weight the ambient term carries
  (c·G_ext = c·4σT_amb⁴ ≡ 4c·σT_amb⁴). It is 2.0 at ε = 1 and 1.64 at
  ε = 0.9; a weight of 1 under-injects the external flux by that factor.
- G_ext = 4σT_amb⁴ is the equilibrium incident radiation from ambient
- G_surf is the P1 incident radiation at the surface
- ε is surface emissivity (1-ε is reflected back into domain)

Derived from the P1 half-range fluxes q_in = G/4 + q/2, q_out = G/4 − q/2
with the surface closure q_in = ε(I_ext + σT_amb⁴) + (1−ε)·q_out (see the
kernel comments above `_marshak_exchange_coefficient`).

This formulation ensures:
1. At equilibrium (G = 4σT⁴, T_amb = T_surf, I_ext = 0): q = 0 ✓
2. External and ambient radiation enter with consistent half-range weights
3. Surface partially reflects (1-ε) of interior radiation back

The production solver applies this relation at the boundary-cell **center**
through the half-cell-corrected coefficient
(`_marshak_effective_coefficient`); this face-value form is the exact
pointwise relation.

# Arguments
- `bc`: Marshak boundary condition
- `G_surf`: Incident radiation at the surface (face) [W/m²]
- `T_surf`: Surface temperature [K] (used only to evaluate ε(T))
- `α_surf`: Surface absorption coefficient [m⁻¹] (not used)
- `ε_surf`: Surface emissivity [-] (not used — the BC's own emissivity governs)
- `t`: Current time [s] (supports dual numbers for AD)

# Returns
- Radiation flux [W/m²], positive = into domain

# Notes
- Type-generic to support automatic differentiation with Dual numbers
- Uses BC emissivity, not ε_surf from prop_cache, for consistency
"""
@inline function radiation_bc_flux(
    bc::MarshakRadiationBC,
    G_surf,
    T_surf,
    α_surf,
    ε_surf,
    t,
)
    c = _marshak_exchange_coefficient(bc, T_surf)
    return c * (_marshak_driving_field(bc, t) - G_surf)
end

"""
	radiation_bc_flux(::ReflectiveRadiationBC, args...)

Zero flux for reflective boundary (∂G/∂n = 0).
"""
@inline function radiation_bc_flux(
    ::ReflectiveRadiationBC,
    G_surf,
    T_surf,
    α_surf,
    ε_surf,
    t,
)
    return zero(t)  # Preserve type for AD
end

"""
	radiation_bc_flux(bc::DirichletRadiationBC, args...)

For Dirichlet BCs, the flux is determined by the interior gradient.
This returns zero as a placeholder - actual implementation needs the interior G value.
"""
@inline function radiation_bc_flux(
    bc::DirichletRadiationBC,
    G_surf,
    T_surf,
    α_surf,
    ε_surf,
    t,
)
    # Dirichlet BC sets G directly; flux determined by interior gradient
    # This is handled specially in the discretization
    return zero(t)  # Preserve type for AD
end

# ============================================================================
# Optical Thickness Diagnostics
# ============================================================================

"""
	optical_thickness(mesh::Mesh1D{NC}, prop_cache::PropertyCache) where NC -> Float64

Compute the optical thickness τ = α × L for the domain.

This is a key parameter for assessing P1 validity:
- τ > 10: Optically thick (P1 excellent)
- 1 < τ < 10: Intermediate (P1 good)  
- 0.1 < τ < 1: Optically thin (P1 fair)
- τ < 0.1: Transparent (P1 poor)

# Arguments
- `mesh`: Mesh with geometry
- `prop_cache`: Property cache with α values

# Returns
- τ: Optical thickness (dimensionless)

# Example
```julia
τ = optical_thickness(mesh, prop_cache)
if τ < 0.1
	@warn "Material is nearly transparent (τ = \$τ). P1 may be inaccurate."
end
```
"""
function optical_thickness(mesh, prop_cache::PropertyCache)
    n_cells = length(mesh.cells)

    # Compute total optical depth: τ = ∫ α dx
    τ_total = 0.0
    @inbounds for i = 1:n_cells
        α = prop_cache.α[i]
        Δx = cell_size(mesh, i)
        τ_total += α * Δx
    end

    return τ_total
end

"""
	optical_thickness_cell(mesh::Mesh1D{NC}, prop_cache::PropertyCache, i::Int) where NC -> Float64

Compute the optical thickness of a single cell: τ_i = α_i × Δx_i

This is useful for identifying which cells are optically thick vs thin.
"""
@inline function optical_thickness_cell(mesh, prop_cache::PropertyCache, i::Int)
    α = prop_cache.α[i]
    Δx = cell_size(mesh, i)
    return α * Δx
end

"""
	check_p1_validity(mesh, prop_cache::PropertyCache; 
					  warn::Bool=true) -> NamedTuple

Check if P1 approximation is valid for the current material state.

Returns a NamedTuple with:
- `τ_total`: Total optical thickness
- `τ_min`: Minimum cell optical thickness
- `τ_max`: Maximum cell optical thickness
- `regime`: String describing optical regime
- `valid`: Bool indicating if P1 is appropriate

# Arguments
- `mesh`: Mesh with geometry
- `prop_cache`: Property cache with α values
- `warn`: If true, print warning for marginal cases

# Example
```julia
validity = check_p1_validity(mesh, prop_cache)
println("Optical regime: \$(validity.regime)")
println("P1 valid: \$(validity.valid)")
```
"""
function check_p1_validity(mesh, prop_cache::PropertyCache; warn::Bool = true)
    n_cells = length(mesh.cells)

    τ_total = optical_thickness(mesh, prop_cache)

    τ_min = Inf
    τ_max = 0.0
    @inbounds for i = 1:n_cells
        τ_i = optical_thickness_cell(mesh, prop_cache, i)
        τ_min = min(τ_min, τ_i)
        τ_max = max(τ_max, τ_i)
    end

    # Classify regime
    if τ_total >= 10
        regime = "optically thick (τ = $(round(τ_total, sigdigits=3)))"
        valid = true
    elseif τ_total > 1
        regime = "intermediate (τ = $(round(τ_total, sigdigits=3)))"
        valid = true
    elseif τ_total > 0.1
        regime = "optically thin (τ = $(round(τ_total, sigdigits=3)))"
        valid = true  # P1 still usable but less accurate
        if warn
            @warn "Material is optically thin (τ = $τ_total). P1 approximation may have reduced accuracy (15-30% error)."
        end
    else
        regime = "transparent (τ = $(round(τ_total, sigdigits=3)))"
        valid = false
        if warn
            @warn "Material is nearly transparent (τ = $τ_total). P1 approximation may be inaccurate (>30% error). Consider discrete ordinates method."
        end
    end

    return (τ_total = τ_total, τ_min = τ_min, τ_max = τ_max, regime = regime, valid = valid)
end

# ============================================================================
# Energy Conservation Diagnostics
# ============================================================================

"""
	radiation_energy_balance(rad_cache::RadiationCache, mesh) -> NamedTuple

Compute radiation energy balance for the domain.

Returns:
- `total_emission`: ∫ Q_emission dV [W/m²] (integrated over depth)
- `total_absorption`: ∫ Q_absorption dV [W/m²]
- `net_source`: total_absorption - total_emission [W/m²]
- `flux_top`: Radiation flux at top boundary [W/m²]
- `flux_bottom`: Radiation flux at bottom boundary [W/m²]
- `balance_error`: (flux_top − flux_bottom) + net_source — zero at a
  converged quasi-steady state

The quasi-steady G equation integrates to ∫(Q_abs − Q_em)dV + ∮q·n̂ dA = 0
with `q_rad` storing the physical +z flux at every face (boundary faces
included), so the discrete residual is `(flux_top − flux_bottom) +
net_source`. The two terms cancel — they do not equal each other.
"""
function radiation_energy_balance(rad_cache::RadiationCache, mesh)
    n_cells = length(mesh.cells)

    total_emission = 0.0
    total_absorption = 0.0

    @inbounds for i = 1:n_cells
        Δx = cell_size(mesh, i)
        total_emission += rad_cache.Q_emission[i] * Δx
        total_absorption += rad_cache.Q_absorption[i] * Δx
    end

    net_source = total_absorption - total_emission

    # Boundary fluxes (top = last face, bottom = first face)
    flux_top = rad_cache.q_rad[end]
    flux_bottom = rad_cache.q_rad[1]

    # Integral of the quasi-steady residual: net_source + ∮q·n̂ = 0
    balance_error = (flux_top - flux_bottom) + net_source

    return (
        total_emission = total_emission,
        total_absorption = total_absorption,
        net_source = net_source,
        flux_top = flux_top,
        flux_bottom = flux_bottom,
        balance_error = balance_error,
    )
end

# P1 validity is a property of the DOMAIN optical thickness τ = ∫α dz — the
# `check_p1_validity(mesh, prop_cache)` method above is the only validity
# check. A per-cell τ_i ≥ 1 criterion is mesh-dependent (refining the mesh
# shrinks every τ_i and degrades the verdict with no change in physics) and
# was removed.
