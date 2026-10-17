# ============================================================================
# Per-cell reference data and tolerance vectors
# ============================================================================

# Smooth floor on the volume used for m/V, as a fraction of the initial
# volume. A cell deactivates long before its volume approaches it.
const VOLUME_FLOOR_FRACTION = 1e-8

# Mass scale of a cell that starts with no condensed mass: its total initial
# mass, or, for an empty cell, 1 kg/m³ over its volume.
@inline function _mass_tolerance_scale(m_solid, m_total, V)
    m_solid > 0 && return m_solid
    m_total > 0 && return m_total
    return 1.0 * V
end

"""
	CellReference(mesh, m0, material; reltol, threshold)

The per-cell reference data of a mesh whose cells hold the initial component
masses `m0` (`NC × n`): initial volumes, the initial condensed mass, the
positive-part scale `δᵢ = reltol · threshold · mˢᵢ(0)` (the cell's mass
absolute tolerance), the volume floor and the heat-capacity floor
`RHO_CP_FLOOR · Vᵢ(0)`.
"""
function CellReference(
    mesh::Mesh{1},
    m0::AbstractMatrix{Float64},
    material::Material{NC};
    reltol::Float64,
    threshold::Float64,
) where {NC}
    topo = mesh.topology
    n = n_cells(topo)
    size(m0) == (NC, n) || throw(DimensionMismatch("m0 must be $NC × $n"))
    x = [p[1] for p in mesh.positions]
    g = MeshGeometry{Float64}(topo)
    geometry!(g, x, topo.cross_section, topo)
    V0 = copy(g.cell_volume)
    all(>(0), V0) || throw(ArgumentError("every initial cell volume must be positive"))
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    m_solid = [sum(m0[j, i] for j = 1:NC if condensed[j]; init = 0.0) for i = 1:n]
    m_total = [sum(m0[:, i]) for i = 1:n]
    scale = [_mass_tolerance_scale(m_solid[i], m_total[i], V0[i]) for i = 1:n]
    δ = reltol * threshold .* scale
    return CellReference(
        V0,
        scale,
        δ,
        VOLUME_FLOOR_FRACTION .* V0,
        RHO_CP_FLOOR .* V0,
        x,
        topo.cross_section,
    )
end

"""
	abstol_vector(layout, setup, m0, T0; reltol, abstol, threshold) -> Vector{Float64}

Absolute tolerance per row. Mass rows `reltol · threshold · mˢᵢ(0)`, with
`mˢᵢ(0)` the cell's total initial condensed mass (never a per-species
inventory, which is zero for a product); temperature rows `abstol` [K];
enthalpy rows `threshold · mˢᵢ(0) · c_p,ref · abstol`, so a cell at the
deactivation threshold resolves its temperature to `abstol` as the
temperature form does.
"""
function abstol_vector(
    layout::StateLayout{NC},
    setup::StateSetup,
    m0::AbstractMatrix{Float64},
    T0::AbstractVector{Float64};
    reltol::Float64,
    abstol::Float64,
    threshold::Float64,
) where {NC}
    tol = fill(abstol, state_length(layout))
    ref = setup.reference
    for i = 1:n_cells(layout)
        for j = 1:NC
            tol[mass_index(layout, i, j)] = ref.solid_mass[i] * reltol * threshold
        end
        tol[energy_index(layout, i)] = _energy_tolerance(
            energy_state(layout),
            setup.laws,
            ref,
            m0,
            T0,
            i,
            abstol,
            threshold,
        )
    end
    return tol
end

@inline _energy_tolerance(::TemperatureState, laws, ref, m0, T0, i, abstol, threshold) =
    abstol

function _energy_tolerance(::EnthalpyState, laws, ref, m0, T0, i, abstol, threshold)
    c = component_heat_capacities(laws, T0[i])
    mtot = sum(view(m0, :, i))
    cp_ref =
        mtot > 0 ? sum(m0[j, i] * c[j] for j in eachindex(c)) / mtot : sum(c) / length(c)
    return threshold * ref.solid_mass[i] * cp_ref * abstol
end
