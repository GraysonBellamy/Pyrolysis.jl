# ============================================================================
# Problem definition and workspace of the finite-volume core
# ============================================================================

"""
	NoSources

No external sources: the residual is the physics alone.
"""
struct NoSources end

"""
	VolumetricSources{FM,FE}

External sources for verification by manufactured solutions, pure functions
of position and time so the Jacobian stays exact:

- `mass(z, t) -> NTuple{NC}`: volumetric mass source of every component
  [kg/(m³·s)], added to the mass rows as `S V`;
- `energy(z, t)`: the temperature-rate source [K/s], added to the energy row as
  `S_T` in the temperature form and as `(Σⱼ mⱼ c_p,ⱼ) S_T` in the inventory form.

Either may be `nothing`. Positions are the cell centers.
"""
struct VolumetricSources{FM,FE}
    mass::FM
    energy::FE
end

"""
	CoreDefinition

The immutable definition the residual is a pure function of, together with
the state vector, the time and the workspace: the state layout, the state
setup (topology, material, reference data, caloric extension), the kernels
in assembler order, the mesh closure (`nothing` on a fixed mesh), and the
external sources.
"""
struct CoreDefinition{L<:StateLayout,S<:StateSetup,K<:Tuple,Cl,Src}
    layout::L
    setup::S
    kernels::K
    closure::Cl
    sources::Src
end

"""
	CoreWorkspace

The mutable buffers of one residual evaluation: the derived state and the
kernel contributions with `Float64` elements, and their twins for every
dual-number element type automatic differentiation evaluates the residual at,
built on first use.
"""
struct CoreWorkspace{DS<:DerivedState,C<:Contributions}
    derived::DS
    contributions::C
    dual::IdDict{DataType,Any}
    T_warm::Vector{Float64}
end

@inline _n_reactions(def::CoreDefinition) = n_reactions(def.setup.material)

"""
	CoreWorkspace(def, T0) -> CoreWorkspace

Buffers for `def`, the inventory inversion warm-started at the cell
temperatures `T0`.
"""
function CoreWorkspace(def::CoreDefinition, T0::AbstractVector{<:Real})
    layout = def.layout
    T_warm = Float64.(T0)
    ds = DerivedState{Float64}(layout, def.setup, T_warm)
    c = Contributions{Float64}(
        n_components(layout),
        _n_reactions(def),
        n_cells(layout),
        n_faces(def.setup.topology),
        n_nodes(layout),
    )
    return CoreWorkspace(ds, c, IdDict{DataType,Any}(), T_warm)
end

"""
	buffers(ws, def, T) -> (derived, contributions)

The derived state and contributions with element type `T`.
"""
@inline buffers(ws::CoreWorkspace, def::CoreDefinition, ::Type{Float64}) =
    (ws.derived, ws.contributions)
function buffers(ws::CoreWorkspace, def::CoreDefinition, ::Type{T}) where {T}
    b = get(ws.dual, T, nothing)
    b === nothing || return b::Tuple{DerivedState,Contributions}
    layout = def.layout
    ds = DerivedState{T}(layout, def.setup, ws.T_warm)
    c = Contributions{T}(
        n_components(layout),
        _n_reactions(def),
        n_cells(layout),
        n_faces(def.setup.topology),
        n_nodes(layout),
    )
    ws.dual[T] = (ds, c)
    return (ds, c)
end
