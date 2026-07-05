"""
Property caching for performance.

Pre-allocated storage for evaluated properties to avoid allocations in inner loops.
"""

"""
	DarcyCache{Tu}

Per-cell and per-face buffers used by Darcy-Fick gas transport. Lives as
`PropertyCache.darcy` and is allocated only when the simulation mode
requires pressure-driven flow; non-Darcy modes carry `nothing` instead so
the per-cell × 4 + per-face × 2 storage is omitted entirely.

# Fields
- `P::Vector{Tu}`     : Pressure field [Pa] (per cell)
- `φ::Vector{Tu}`     : Porosity [-] (per cell)
- `κ::Vector{Tu}`     : Permeability [m²] (per cell)
- `μ::Vector{Tu}`     : Gas viscosity [Pa·s] (per cell)
- `v_gas::Vector{Tu}` : Darcy velocity [m/s] (per face)
- `∇P::Vector{Tu}`    : Pressure gradient [Pa/m] (per face)
"""
mutable struct DarcyCache{Tu}
    P::Vector{Tu}
    φ::Vector{Tu}
    κ::Vector{Tu}
    μ::Vector{Tu}
    v_gas::Vector{Tu}
    ∇P::Vector{Tu}
end

DarcyCache{Tu}(n_cells::Int, n_faces::Int) where {Tu} = DarcyCache{Tu}(
    zeros(Tu, n_cells),
    zeros(Tu, n_cells),
    zeros(Tu, n_cells),
    zeros(Tu, n_cells),
    zeros(Tu, n_faces),
    zeros(Tu, n_faces),
)

function resize_darcy_cache!(c::DarcyCache, n_cells::Int, n_faces::Int)
    resize!(c.P, n_cells)
    resize!(c.φ, n_cells)
    resize!(c.κ, n_cells)
    resize!(c.μ, n_cells)
    resize!(c.v_gas, n_faces)
    resize!(c.∇P, n_faces)
    return c
end

"""
	PropertyCache{Tu, DarcyT}

Pre-allocated storage for evaluated material properties. Parameterized
on the live scalar type `Tu` so the AD workspace can hold
`ForwardDiff.Dual`-typed property values without re-allocation. The
`DarcyT` parameter is either `DarcyCache{Tu}` (when pressure-driven gas
transport is enabled) or `Nothing` (Fickian-only modes), so the Darcy
buffers are allocated only when actually needed.

# Fields
## Per-cell properties
- `ρ::Vector{Tu}`: Mixture density [kg/m³]
- `ρcp::Vector{Tu}`: Volumetric heat capacity [J/(m³·K)]
- `k::Vector{Tu}`: Thermal conductivity [W/(m·K)]
- `λ::Vector{Tu}`: Gas transfer coefficient [m²/s]
- `ε::Vector{Tu}`: Emissivity [-]
- `α::Vector{Tu}`: Absorption coefficient [m⁻¹]

## Per-face properties
- `k_face::Vector{Tu}`: Face-centered conductivity [W/(m·K)]
- `λ_face::Vector{Tu}`: Face-centered gas transfer coefficient [m²/s]

## Optional Darcy sub-cache
- `darcy::DarcyT`: `DarcyCache{Tu}` for pressure-driven transport, or
  `Nothing` for Fickian-only configurations.
"""
mutable struct PropertyCache{Tu,DarcyT}
    # Per-cell properties
    ρ::Vector{Tu}
    ρcp::Vector{Tu}
    k::Vector{Tu}
    λ::Vector{Tu}
    ε::Vector{Tu}
    α::Vector{Tu}

    # Per-face properties (for flux computation)
    k_face::Vector{Tu}
    λ_face::Vector{Tu}

    # Optional Darcy sub-cache
    darcy::DarcyT
end


"""
	PropertyCache{Tu}(n_cells::Int, n_faces::Int; with_darcy::Bool=true) where {Tu}
	PropertyCache(n_cells::Int, n_faces::Int)         # Tu = Float64

Create a PropertyCache with pre-allocated arrays. Pass `with_darcy=false`
to skip the Darcy sub-cache for Fickian-only configurations.
"""
function PropertyCache{Tu}(n_cells::Int, n_faces::Int; with_darcy::Bool = true) where {Tu}
    # Branch on `with_darcy` here (not via `typeof(darcy)` on a Union value) so the
    # function barrier below sees a concrete `darcy` type in each path — keeps the
    # inner constructor call free of a `Union{DarcyCache,Nothing}` (which JET reads
    # as a possible `convert(::DarcyCache, ::Nothing)`).
    if with_darcy
        return _build_property_cache(Tu, n_cells, n_faces, DarcyCache{Tu}(n_cells, n_faces))
    else
        return _build_property_cache(Tu, n_cells, n_faces, nothing)
    end
end

function _build_property_cache(
    ::Type{Tu},
    n_cells::Int,
    n_faces::Int,
    darcy::D,
) where {Tu,D}
    return PropertyCache{Tu,D}(
        zeros(Tu, n_cells),  # ρ
        zeros(Tu, n_cells),  # ρcp
        zeros(Tu, n_cells),  # k
        zeros(Tu, n_cells),  # λ
        zeros(Tu, n_cells),  # ε
        zeros(Tu, n_cells),  # α
        zeros(Tu, n_faces),  # k_face
        zeros(Tu, n_faces),  # λ_face
        darcy,
    )
end

PropertyCache(n_cells::Int, n_faces::Int) = PropertyCache{Float64}(n_cells, n_faces)

"""
	PropertyCache(mesh::Mesh1D)

Create a Float64-typed PropertyCache sized for the given mesh.
"""
PropertyCache(mesh) = PropertyCache(length(mesh.cells), length(mesh.faces))

"""
	resize_cache!(cache::PropertyCache, n_cells::Int, n_faces::Int)

Resize cache arrays (e.g., after mesh refinement). The Darcy sub-cache,
if present, is resized through the same call.
"""
function resize_cache!(cache::PropertyCache, n_cells::Int, n_faces::Int)
    resize!(cache.ρ, n_cells)
    resize!(cache.ρcp, n_cells)
    resize!(cache.k, n_cells)
    resize!(cache.λ, n_cells)
    resize!(cache.ε, n_cells)
    resize!(cache.α, n_cells)
    resize!(cache.k_face, n_faces)
    resize!(cache.λ_face, n_faces)
    if cache.darcy !== nothing
        resize_darcy_cache!(cache.darcy, n_cells, n_faces)
    end
    return cache
end

"""
	update_properties!(cache::PropertyCache, state, mesh, material::Material{NC}) where NC

Evaluate all material properties at the current state and store in cache.

State is read from `state::StateCache` (T, ξ, cell-center / face positions);
`mesh` is consulted only for face connectivity (`mesh.faces[f].cell_ids`).
The function is allocation-free after the cache is sized to match.
"""
function update_properties!(
    cache::PropertyCache,
    state,
    mesh,
    material::Material{NC},
) where {NC}
    n_cells_val = length(state.T)
    n_faces_val = length(state.A_face)

    # Ensure cache is properly sized
    if length(cache.ρ) != n_cells_val
        resize_cache!(cache, n_cells_val, n_faces_val)
    end

    # Evaluate cell-centered properties.
    # volume_fractions / solid_volume_fractions are computed
    # once per cell and threaded into the effective-property helpers, instead
    # of being recomputed inside each call.
    # The previous step's `cache.k[i]` warm-starts the
    # Bruggeman Newton iteration when the conductivity mixing rule is BRUGGEMAN.
    @inbounds for i = 1:n_cells_val
        T_i = state.T[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))

        v_full = volume_fractions(ξ_i, T_i, material)
        # Depleted-cell floor: falls back to v_full when no
        # matrix volume remains, so k_eff degrades to the gas-mixture value.
        v_solid = mixing_volume_fractions(ξ_i, T_i, material)
        k_prev = cache.k[i]

        cache.ρ[i] = mixture_density(ξ_i, T_i, material)
        cache.ρcp[i] = effective_heat_capacity(ξ_i, T_i, material)
        cache.k[i] = _conductivity_from_v(v_solid, T_i, material; k_init = k_prev)
        cache.λ[i] = _gas_transfer_from_v(v_full, T_i, material)
        cache.ε[i] = _emissivity_from_v(v_solid, material)
        cache.α[i] = effective_absorption(ξ_i, material)
    end

    # Compute face-centered properties using distance-weighted harmonic mean
    # (resistance-in-series for non-uniform meshes).
    @inbounds for f = 1:n_faces_val
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids
        face_pos = state.face_positions[f]

        if left_id == 0
            # Left boundary: use right cell value
            cache.k_face[f] = cache.k[right_id]
            cache.λ_face[f] = cache.λ[right_id]
        elseif right_id == 0
            # Right boundary: use left cell value
            cache.k_face[f] = cache.k[left_id]
            cache.λ_face[f] = cache.λ[left_id]
        else
            d_L = abs(face_pos - state.cell_centers[left_id])
            d_R = abs(state.cell_centers[right_id] - face_pos)

            cache.k_face[f] =
                harmonic_mean_weighted(cache.k[left_id], cache.k[right_id], d_L, d_R)
            cache.λ_face[f] =
                harmonic_mean_weighted(cache.λ[left_id], cache.λ[right_id], d_L, d_R)
        end
    end

    return cache
end
