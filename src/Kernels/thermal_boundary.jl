# ============================================================================
# Thermal boundary conditions
# ============================================================================

"""
	BoundaryPatch{BC}

The boundary faces carrying one tag and the condition applied on them.
"""
struct BoundaryPatch{BC}
    tag::Symbol
    faces::Vector{Int}
    bc::BC
end

"""
	ThermalBoundaryKernel{P<:Tuple} <: Kernel

The energy flux through every boundary face, one [`BoundaryPatch`](@ref) per
boundary tag. The flux is written along the outward face normal.
"""
struct ThermalBoundaryKernel{P<:Tuple} <: Kernel
    patches::P
    d_floor::Vector{Float64}
end

"""
	ThermalBoundaryKernel(topology, geometry, bc_set)

One patch per boundary tag of the mesh, holding the thermal condition the
boundary-condition set assigns to that tag.
"""
function ThermalBoundaryKernel(
    topo::MeshTopology,
    g::MeshGeometry,
    bc_set::BoundaryConditionSet,
)
    tags = boundary_tags(topo)
    for tag in keys(bc_set.thermal)
        tag in tags || throw(
            ArgumentError(
                "the boundary conditions name the tag :$tag, which the mesh does not have " *
                "(mesh tags: $tags)",
            ),
        )
    end
    patches = Tuple(
        BoundaryPatch(tag, boundary_faces(topo, tag), bc_set.thermal[tag]) for tag in tags
    )
    d = zeros(n_faces(topo))
    for f = 1:n_faces(topo)
        o = topo.faces[f].cells[1]
        d[f] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, f, o)
    end
    return ThermalBoundaryKernel(patches, d)
end

function contribute!(
    c::Contributions,
    kernel::ThermalBoundaryKernel,
    ds::DerivedState,
    topo::MeshTopology,
    material,
    t,
)
    _patches!(c, kernel.patches, kernel, ds, topo, t)
    return c
end

# Peel the patch tuple so each patch's condition is a concrete type.
@inline _patches!(c, ::Tuple{}, kernel, ds, topo, t) = nothing
@inline function _patches!(c, patches::Tuple, kernel, ds, topo, t)
    p = first(patches)
    @inbounds for f in p.faces
        c.F_energy[f] += boundary_energy_flux(p.bc, f, kernel, ds, topo, t)
    end
    return _patches!(c, Base.tail(patches), kernel, ds, topo, t)
end

"""
	boundary_energy_flux(bc, f, kernel, ds, topology, t)

Outward energy flux per unit area [W/m²] through boundary face `f` under the
thermal condition `bc`.
"""
@inline boundary_energy_flux(::AdiabaticBC, f, kernel, ds, topo, t) = zero(eltype(ds.T))

@inline boundary_energy_flux(bc::HeatFluxBC, f, kernel, ds, topo, t) =
    -evaluate_param(bc.flux, t) + zero(eltype(ds.T))

@inline function boundary_energy_flux(bc::DirichletThermalBC, f, kernel, ds, topo, t)
    i = topo.faces[f].cells[1]
    d = smooth_floor(center_to_face_distance(ds.geometry, topo, f, i), kernel.d_floor[f])
    return ds.k[i] * (ds.T[i] - evaluate_param(bc.T, t)) / d
end
