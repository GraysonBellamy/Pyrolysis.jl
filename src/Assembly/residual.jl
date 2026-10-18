# ============================================================================
# Residual
# ============================================================================

"""
	residual!(du, u, def::CoreDefinition, t, ws::CoreWorkspace) -> nothing

The right-hand side of the core's mass-matrix system `M du/dt = f(u, t)`: the
state mapping, the kernels in assembler order, the external sources, then the
row assembly. A pure function of `(u, t)` given `def`; `ws` holds buffers
only. Dual-number states evaluate in their own buffers.
"""
function residual!(
    du::AbstractVector,
    u::AbstractVector,
    def::CoreDefinition,
    t,
    ws::CoreWorkspace,
)
    ds, c = buffers(ws, def, eltype(u))
    _residual!(du, u, def, t, ds, c)
    return nothing
end

function _residual!(du, u, def::CoreDefinition, t, ds::DerivedState, c::Contributions)
    setup = def.setup
    derive!(ds, u, def.layout, setup)
    reset!(c)
    _kernels!(c, def.kernels, ds, setup.topology, setup.material, t)
    add_sources!(c, def.sources, ds, setup, t)
    assemble_rows!(du, u, c, ds, def.layout, setup, def.sources, t)
    return nothing
end

@inline _kernels!(c, ::Tuple{}, ds, topo, material, t) = nothing
@inline function _kernels!(c, kernels::Tuple, ds, topo, material, t)
    contribute!(c, first(kernels), ds, topo, material, t)
    return _kernels!(c, Base.tail(kernels), ds, topo, material, t)
end

# ============================================================================
# External sources
# ============================================================================

"""
	add_sources!(c, sources, ds, setup, t)

Add the external sources to the contributions: mass sources to `S_mass`; the
energy source is applied by the row assembly, which knows the energy state.
"""
@inline add_sources!(c, ::NoSources, ds, setup, t) = nothing

function add_sources!(c::Contributions{NC}, s::VolumetricSources, ds, setup, t) where {NC}
    s.mass === nothing && return nothing
    g = ds.geometry
    @inbounds for i in axes(c.S_mass, 2)
        S = s.mass(g.cell_center[i][1], t)
        for j = 1:NC
            c.S_mass[j, i] += S[j] * ds.V_eval[i]
        end
    end
    return nothing
end

# The external temperature-rate source of cell `i` [K/s].
@inline _energy_source(::NoSources, ds, i, t) = zero(eltype(ds.T))
@inline _energy_source(s::VolumetricSources, ds, i, t) =
    s.energy === nothing ? zero(eltype(ds.T)) : s.energy(ds.geometry.cell_center[i][1], t)

# ============================================================================
# Row assembly
# ============================================================================

"""
	assemble_rows!(du, u, c, ds, layout, setup, sources, t)

Write the row derivatives from the kernel contributions: the only writer of
`du`. Every face flux is scattered as `∓F_f A_f` into its owner and neighbour;
the mass rows take their reaction and face-flux parts; the energy row is

- temperature form: `dTᵢ/dt = (Eᵢ − Σⱼ hⱼ(Tᵢ) Ṁᵢⱼ) / Σⱼ mᵢⱼ c_p,ⱼ`,
- inventory form: `dHᵢ/dt = Eᵢ + Σⱼ ωᵢⱼ Vᵢ hⱼ(Tᵢ)`,

with `Eᵢ` the energy sources plus the face energy fluxes and `Ṁᵢⱼ` the
face-flux part of the mass rate. An external temperature-rate source `S_T`
adds `S_T` to the temperature row and `(Σⱼ mᵢⱼ c_p,ⱼ) S_T` to the inventory row.
The face-intensity rows take the residuals of their algebraic equations.
Rows the layout marks dead are identity rows.
"""
function assemble_rows!(
    du::AbstractVector,
    u::AbstractVector,
    c::Contributions{NC},
    ds::DerivedState{NC},
    layout::StateLayout{NC},
    setup::StateSetup,
    sources,
    t,
) where {NC}
    topo = setup.topology
    g = ds.geometry
    n = n_cells(layout)

    @inbounds for i = 1:n
        for j = 1:NC
            c.M_flux[j, i] = zero(eltype(c.M_flux))
        end
        du[energy_index(layout, i)] = c.S_energy[i]
    end
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        A = g.face_area[f]
        q = c.F_energy[f] * A
        du[energy_index(layout, o)] -= q
        nb == 0 || (du[energy_index(layout, nb)] += q)
        for j = 1:NC
            F = c.F_mass[j, f] * A
            c.M_flux[j, o] -= F
            nb == 0 || (c.M_flux[j, nb] += F)
        end
    end
    @inbounds for i = 1:n
        for j = 1:NC
            du[mass_index(layout, i, j)] = c.S_mass[j, i] + c.M_flux[j, i]
        end
        E = du[energy_index(layout, i)]
        S_T = _energy_source(sources, ds, i, t)
        du[energy_index(layout, i)] = _energy_row(energy_state(layout), E, S_T, c, ds, i)
    end
    if has_intensities(layout)
        @inbounds for f = 1:n_faces(topo)
            du[intensity_index(layout, f)] = c.r_intensity[f]
        end
    end
    _dead_rows!(du, u, layout)
    return nothing
end

@inline function _energy_row(
    ::TemperatureState,
    E,
    S_T,
    c::Contributions{NC},
    ds,
    i,
) where {NC}
    for j = 1:NC
        E -= ds.h[j, i] * c.M_flux[j, i]
    end
    return E / ds.cp_sum[i] + S_T
end

@inline function _energy_row(
    ::EnthalpyState,
    E,
    S_T,
    c::Contributions{NC},
    ds,
    i,
) where {NC}
    for j = 1:NC
        E += c.S_mass[j, i] * ds.h[j, i]
    end
    return E + ds.cp_sum[i] * S_T
end

# Dead rows: zero derivative on the differential rows, `u` itself on the
# algebraic rows, so each is an identity row of the Newton matrix.
function _dead_rows!(du, u, layout::StateLayout)
    mask = layout.mask
    all(mask) && return nothing
    diff = layout.differential
    @inbounds for k in eachindex(mask)
        mask[k] && continue
        du[k] = diff[k] ? zero(eltype(du)) : u[k]
    end
    return nothing
end
