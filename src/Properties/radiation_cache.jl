"""
Radiation property caching for P1 radiation transport.

Pre-allocated storage for radiation-specific properties to avoid allocations
in inner loops. The RadiationCache complements the main PropertyCache by
storing the radiation intensity field G and derived quantities.

Reference: P1_Volumetric_Reradiation_Implementation.md Section 4.2
"""

# ============================================================================
# RadiationCache Structure
# ============================================================================

"""
	RadiationCache

Pre-allocated storage for P1 radiation properties and state.

The cache stores both the radiation intensity field G (the primary P1 variable)
and derived quantities needed for flux computation.

# Fields

## Radiation Diffusion Coefficient
- `D_rad::Vector{Float64}`: Cell-centered radiation diffusion coefficient [m].
  D_rad = 1/(3α) where α is absorption coefficient.
  
- `D_rad_face::Vector{Float64}`: Face-centered D_rad [m], computed as 
  distance-weighted harmonic mean for accurate interface treatment.

## Radiation State
- `G::Vector{Float64}`: Incident radiation intensity [W/m²]. This is the
  primary variable solved in the P1 equation.
  
- `q_rad::Vector{Float64}`: Radiative heat flux at faces [W/m²].
  q_rad = -D_rad × ∂G/∂x

## Source Terms
- `Q_emission::Vector{Float64}`: Volumetric emission source [W/m³].
  Q_emission = 4σαT⁴ (blackbody emission into the radiation field)
  
- `Q_absorption::Vector{Float64}`: Volumetric absorption [W/m³].
  Q_absorption = αG (energy removed from radiation field)

# Usage
```julia
# Create cache for 100-cell mesh
rad_cache = RadiationCache(100)

# Or create from mesh
rad_cache = RadiationCache(mesh)

# Initialize G field to ambient radiation
initialize_radiation_intensity!(rad_cache, 300.0)

# Update properties after state change
update_radiation_properties!(rad_cache, mesh, material, prop_cache)
```

# Memory Layout
Arrays are sized for n_cells (cell-centered) and n_faces (face-centered):
- Cell arrays: D_rad, G, Q_emission, Q_absorption (length = n_cells)
- Face arrays: D_rad_face, q_rad (length = n_faces = n_cells + 1 for 1D)

# See Also
- `PropertyCache` for thermal and transport properties
- `update_radiation_properties!` for property evaluation
- `compute_radiation_fluxes!` for flux computation
"""
mutable struct RadiationCache
    # Cell-centered diffusion coefficient
    D_rad::Vector{Float64}       # Radiation diffusion coefficient [m]

    # Face-centered diffusion coefficient
    D_rad_face::Vector{Float64}  # Face-averaged D_rad [m]

    # Radiation intensity field (primary P1 variable)
    G::Vector{Float64}           # Incident radiation intensity [W/m²]

    # Radiation flux at faces
    q_rad::Vector{Float64}       # Radiative heat flux [W/m²]

    # Volumetric source terms
    Q_emission::Vector{Float64}  # Emission source: 4σαT⁴ [W/m³]
    Q_absorption::Vector{Float64} # Absorption: αG [W/m³]

    # Additional fields for P1 integration
    α::Vector{Float64}           # Volumetric absorption coefficient [1/m]
    S_rad::Vector{Float64}       # Net radiation source α(G - 4σT⁴) [W/m³], positive = heating

    # Pre-allocated work arrays for Newton iteration (avoid allocations in hot loop)
    newton_R::Vector{Float64}    # Residual vector
    newton_ΔG::Vector{Float64}   # Update vector
    newton_dl::Vector{Float64}   # Lower diagonal of tridiagonal Jacobian
    newton_d::Vector{Float64}    # Main diagonal of tridiagonal Jacobian
    newton_du::Vector{Float64}   # Upper diagonal of tridiagonal Jacobian
    newton_d_work::Vector{Float64}  # Work array for Thomas algorithm (preserves newton_d)
end

# ============================================================================
# Constructors
# ============================================================================

"""
	RadiationCache(n_cells::Int)

Create a RadiationCache for a mesh with n_cells cells.

For a 1D mesh, the number of faces is n_cells + 1.
"""
function RadiationCache(n_cells::Int)
    n_faces = n_cells + 1
    RadiationCache(
        zeros(n_cells),      # D_rad
        zeros(n_faces),      # D_rad_face
        zeros(n_cells),      # G
        zeros(n_faces),      # q_rad
        zeros(n_cells),      # Q_emission
        zeros(n_cells),      # Q_absorption
        zeros(n_cells),      # α
        zeros(n_cells),      # S_rad
        zeros(n_cells),      # newton_R
        zeros(n_cells),      # newton_ΔG
        zeros(n_cells - 1),  # newton_dl
        zeros(n_cells),      # newton_d
        zeros(n_cells - 1),  # newton_du
        zeros(n_cells),      # newton_d_work
    )
end

"""
	RadiationCache(mesh::Mesh1D)

Create a RadiationCache sized for the given mesh.
"""
RadiationCache(mesh) = RadiationCache(length(mesh.cells))

# ============================================================================
# Resize Function
# ============================================================================

"""
	resize_radiation_cache!(cache::RadiationCache, n_cells::Int)

Resize radiation cache arrays (e.g., after mesh refinement).

This function maintains the allocated arrays but changes their size.
Contents are preserved where possible but may need reinitialization.
"""
function resize_radiation_cache!(cache::RadiationCache, n_cells::Int)
    n_faces = n_cells + 1

    resize!(cache.D_rad, n_cells)
    resize!(cache.D_rad_face, n_faces)
    resize!(cache.G, n_cells)
    resize!(cache.q_rad, n_faces)
    resize!(cache.Q_emission, n_cells)
    resize!(cache.Q_absorption, n_cells)
    resize!(cache.α, n_cells)
    resize!(cache.S_rad, n_cells)
    resize!(cache.newton_R, n_cells)
    resize!(cache.newton_ΔG, n_cells)
    resize!(cache.newton_dl, n_cells - 1)
    resize!(cache.newton_d, n_cells)
    resize!(cache.newton_du, n_cells - 1)
    resize!(cache.newton_d_work, n_cells)

    return cache
end

# ============================================================================
# Initialization Functions
# ============================================================================

"""
	initialize_radiation_intensity!(cache::RadiationCache, T_ambient::Float64)

Initialize the radiation intensity field G to equilibrium with ambient temperature.

For a medium in thermal equilibrium at temperature T, the incident radiation is:
	G = 4σT⁴ (blackbody radiation field)

In the standard P1 approximation, the equilibrium G is the blackbody value.
Emissivity only appears at boundaries, not in the volumetric source term.

# Arguments
- `cache`: RadiationCache to initialize
- `T_ambient`: Ambient/initial temperature [K]
"""
function initialize_radiation_intensity!(cache::RadiationCache, T_ambient::Float64)
    # Equilibrium radiation field: G = 4σT⁴ (blackbody)
    G_ambient = 4.0 * STEFAN_BOLTZMANN * T_ambient^4
    fill!(cache.G, G_ambient)

    # Zero out fluxes and sources (will be computed on first update)
    fill!(cache.q_rad, 0.0)
    fill!(cache.Q_emission, 0.0)
    fill!(cache.Q_absorption, 0.0)

    return cache
end

"""
	initialize_radiation_intensity!(cache::RadiationCache, mesh::Mesh1D)

Initialize G field based on cell temperatures in the mesh.

Uses the local cell temperature to set initial G = 4σT⁴ for each cell.
"""
function initialize_radiation_intensity!(cache::RadiationCache, mesh)
    n_cells = length(mesh.cells)

    @inbounds for i = 1:n_cells
        T = mesh.cell_states[i].T
        cache.G[i] = 4.0 * STEFAN_BOLTZMANN * T^4
    end

    fill!(cache.q_rad, 0.0)
    fill!(cache.Q_emission, 0.0)
    fill!(cache.Q_absorption, 0.0)

    return cache
end

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
	radiation_energy_density(cache::RadiationCache, i::Int) -> Float64

Get the radiation energy density at cell i.

In the P1 approximation, the radiation energy density is:
	E_rad = G / c

where c is the speed of light. This is typically very small compared to
thermal energy but useful for diagnostics.
"""
@inline function radiation_energy_density(cache::RadiationCache, i::Int)
    c_light = 2.998e8  # m/s
    return cache.G[i] / c_light
end

"""
	net_radiation_source(cache::RadiationCache, i::Int) -> Float64

Get the net volumetric radiation source for cell i.

Net source = Q_absorption - Q_emission = αG - 4σαT⁴

Positive means net absorption (heating), negative means net emission (cooling).
"""
@inline function net_radiation_source(cache::RadiationCache, i::Int)
    return cache.Q_absorption[i] - cache.Q_emission[i]
end

"""
	total_radiation_flux(cache::RadiationCache) -> Float64

Get the total radiative flux through the domain.

Returns the sum of absolute radiation fluxes, useful for energy balance checks.
"""
function total_radiation_flux(cache::RadiationCache)
    return sum(abs, cache.q_rad)
end
