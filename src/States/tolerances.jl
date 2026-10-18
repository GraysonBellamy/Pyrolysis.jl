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
	abstol_vector(layout, setup, m0, T0; reltol, abstol, threshold,
	              intensity_scale = 1.0, time_scale = 1.0) -> Vector{Float64}

Absolute tolerance per row. Mass rows `reltol · threshold · mˢᵢ(0)`, with
`mˢᵢ(0)` the cell's total initial condensed mass (never a per-species
inventory, which is zero for a product); temperature rows `abstol` [K];
enthalpy rows `threshold · mˢᵢ(0) · c_p,ref · abstol`, so a cell at the
deactivation threshold resolves its temperature to `abstol` as the
temperature form does; face-intensity rows `reltol · intensity_scale`, the
scale being the largest incident intensity [W/m²].

On a moving mesh: node positions `reltol · threshold · hₙ`, with `hₙ` the
thinnest initial cell beside node `n`, so a position error moves a cell's
concentrations by the relative amount its mass tolerance allows; node
velocities that over `time_scale` [s], the first step; the history pair the
mass tolerance of its cell; the cross-section measure `reltol · A₀`; and its
lateral rate the velocity tolerance spread over the initial thickness.
"""
function abstol_vector(
    layout::StateLayout{NC},
    setup::StateSetup,
    m0::AbstractMatrix{Float64},
    T0::AbstractVector{Float64};
    reltol::Float64,
    abstol::Float64,
    threshold::Float64,
    intensity_scale::Float64 = 1.0,
    time_scale::Float64 = 1.0,
) where {NC}
    tol = fill(abstol, state_length(layout))
    if has_intensities(layout)
        for f = 1:n_nodes(layout)
            tol[intensity_index(layout, f)] = reltol * intensity_scale
        end
    end
    ref = setup.reference
    if is_moving(layout)
        x = ref.positions
        nn = n_nodes(layout)
        w_min = Inf
        for nd = 1:nn
            h = min(nd > 1 ? x[nd] - x[nd-1] : Inf, nd < nn ? x[nd+1] - x[nd] : Inf)
            tol[position_index(layout, nd)] = reltol * threshold * h
            tol[velocity_index(layout, nd)] = reltol * threshold * h / time_scale
            w_min = min(w_min, tol[velocity_index(layout, nd)])
        end
        if has_chi(layout)
            tol[area_index(layout)] = reltol * ref.cross_section
            tol[lateral_rate_index(layout)] = w_min / (x[end] - x[1])
        end
    end
    for i = 1:n_cells(layout)
        for j = 1:NC
            tol[mass_index(layout, i, j)] = ref.solid_mass[i] * reltol * threshold
        end
        if has_chi(layout)
            for k = 1:2
                tol[history_index(layout, i, k)] = ref.solid_mass[i] * reltol * threshold
            end
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
