# ============================================================================
# Residual
# ============================================================================

"""
	residual!(du, u, def::CoreDefinition, t, ws::CoreWorkspace) -> nothing

The right-hand side of the core's mass-matrix system `M du/dt = f(u, t)`: the
state mapping, the node corrections of the mesh closure, the kernels in
assembler order, the mesh closure (strain rates from the reaction sources,
velocity rows), the external sources, then the row assembly. A pure function of `(u, t)` given `def`; `ws` holds buffers
only. Dual-number states evaluate in their own buffers, chosen by the state's
element type, so `t` must be a plain number: a forward-mode derivative in
time alone (a dual `t` with a `Float64` state) is not supported, and every
integrator policy takes its time derivative by finite differences.
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
    correction!(ds, def.closure, setup)
    reset!(c)
    _kernels!(c, def.kernels, ds, setup.topology, setup.material, t)
    closure_rows!(c, def.closure, ds, setup.topology, setup.material, t)
    add_sources!(c, def.sources, ds, setup, t)
    assemble_rows!(du, u, c, ds, def.layout, setup, def.closure, def.sources, t)
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
	assemble_rows!(du, u, c, ds, layout, setup, closure, sources, t)

Write the row derivatives from the kernel contributions: the only writer of
`du`. Every face flux is scattered as `∓F_f A_f` into its owner and neighbour;
the mass rows take their reaction and face-flux parts; the energy row is

- temperature form: `dTᵢ/dt = (Eᵢ − Σⱼ hⱼ(Tᵢ) Ṁᵢⱼ) / Σⱼ mᵢⱼ c_p,ⱼ`,
- inventory form: `dHᵢ/dt = Eᵢ + Σⱼ ωᵢⱼ Vᵢ hⱼ(Tᵢ)`,

with `Eᵢ` the energy sources plus the face energy fluxes and `Ṁᵢⱼ` the
face-flux part of the mass rate. An external temperature-rate source `S_T`
adds `S_T` to the temperature row and `(Σⱼ mᵢⱼ c_p,ⱼ) S_T` to the inventory row.
The face-intensity rows take the residuals of their algebraic equations.

On a moving mesh the position rows are `dxₙ/dt = wₙ` and the velocity rows
take the closure's residuals; under a variable cross section the history
rows take their sources and face fluxes, and the area and lateral-rate rows
the residuals of [`lateral_rows`](@ref) at the column sums of the history
pair and its rates. Rows the layout marks dead are identity rows.
"""
function assemble_rows!(
    du::AbstractVector,
    u::AbstractVector,
    c::Contributions{NC},
    ds::DerivedState{NC},
    layout::StateLayout{NC},
    setup::StateSetup,
    closure,
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
    if is_moving(layout)
        @inbounds for nd = 1:n_nodes(layout)
            du[position_index(layout, nd)] = ds.w[nd]
            du[velocity_index(layout, nd)] = c.r_velocity[nd]
        end
    end
    if has_chi(layout)
        _history_rows!(du, c, ds, layout, topo)
        _lateral_rows!(du, ds, layout, closure.lateral)
    end
    _dead_rows!(du, u, layout)
    return nothing
end

# The history rows: each cell's source plus the scattered face fluxes.
function _history_rows!(du, c::Contributions, ds, layout::StateLayout, topo)
    g = ds.geometry
    @inbounds for i = 1:n_cells(layout), k = 1:2
        du[history_index(layout, i, k)] = c.S_carried[k, i]
    end
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        A = g.face_area[f]
        for k = 1:2
            F = c.F_carried[k, f] * A
            du[history_index(layout, o, k)] -= F
            nb == 0 || (du[history_index(layout, nb, k)] += F)
        end
    end
    return nothing
end

# The area and lateral-rate rows at the column sums of the live cells'
# history pair and its rates.
function _lateral_rows!(du, ds, layout::StateLayout, lat::LateralShrinkage)
    T = eltype(du)
    ΣG, ΣM, ΣĠ, ΣṀ = zero(T), zero(T), zero(T), zero(T)
    @inbounds for i = 1:n_cells(layout)
        kG, kM = history_index(layout, i, 1), history_index(layout, i, 2)
        layout.mask[kG] || continue
        ΣG += ds.G[i]
        ΣM += ds.M[i]
        ΣĠ += du[kG]
        ΣṀ += du[kM]
    end
    r = lateral_rows(lat, SVector(ΣG, ΣM, ΣĠ, ΣṀ, ds.θ_A[], ds.A[]))
    du[area_index(layout)] = r[1]
    du[lateral_rate_index(layout)] = r[2]
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
