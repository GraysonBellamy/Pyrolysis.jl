# ============================================================================
# Kernel contributions
# ============================================================================

"""
	Kernel

A physics kernel. [`contribute!`](@ref) adds its face fluxes and cell sources
to a [`Contributions`](@ref) buffer, reading only the derived state, the mesh
topology and metrics, the material and the boundary conditions it holds; it
never reads the state vector, the row derivatives or another kernel's output
except through the buffer, in the order the assembler fixes.
"""
abstract type Kernel end

"""
	Contributions{NC,T}

What the kernels return, with the derived state's element type `T`:

- `F_mass`: mass flux of every component per unit face area along the face
  normal (`NC × n_faces`) [kg/(m²·s)]
- `F_energy`: energy flux per unit face area along the face normal — conduction
  and the enthalpy carried by the mass fluxes (`n_faces`) [W/m²]
- `S_mass`: reaction mass source of every component, `ωᵢⱼ Vᵢ` (`NC × n_cells`) [kg/s]
- `S_energy`: cell energy source — reaction heat and absorbed radiation, without
  the formation term (`n_cells`) [W]
- `θ`: cell strain rate (`n_cells`) [1/s]
- `M_flux`: the face-flux part of every mass rate, `−Σ_f σ_if F_f A_f`
  (`NC × n_cells`) [kg/s], filled by the assembler
- `T_face`: face temperatures (`n_faces`) [K]: on a boundary face the surface
  temperature of its energy balance
- `r_intensity`: the residual of the algebraic row of every face intensity
  (`n_faces`) [W/m²]
- `r_velocity`: the residual of the algebraic row of every node velocity
  (`n_nodes`), written by the mesh closure
- `F_carried`: flux of the conversion-history pair `G`, `M` per unit face area
  along the face normal (`2 × n_faces`) [kg/(m²·s)]
- `S_carried`: cell source of the history pair (`2 × n_cells`) [kg/s]

`rates`, `heats` and `species` are per-cell scratch of the reaction kernel.
"""
struct Contributions{NC,T}
    F_mass::Matrix{T}
    F_energy::Vector{T}
    S_mass::Matrix{T}
    S_energy::Vector{T}
    θ::Vector{T}
    M_flux::Matrix{T}
    T_face::Vector{T}
    r_intensity::Vector{T}
    r_velocity::Vector{T}
    F_carried::Matrix{T}
    S_carried::Matrix{T}
    rates::Vector{T}
    heats::Vector{T}
    species::Vector{T}
end

"""
	Contributions{T}(NC, NR, n_cells, n_faces, n_nodes)

Zero-filled contribution buffers.
"""
function Contributions{T}(
    NC::Int,
    NR::Int,
    n_cells::Int,
    n_faces::Int,
    n_nodes::Int,
) where {T}
    return Contributions{NC,T}(
        zeros(T, NC, n_faces),
        zeros(T, n_faces),
        zeros(T, NC, n_cells),
        zeros(T, n_cells),
        zeros(T, n_cells),
        zeros(T, NC, n_cells),
        zeros(T, n_faces),
        zeros(T, n_faces),
        zeros(T, n_nodes),
        zeros(T, 2, n_faces),
        zeros(T, 2, n_cells),
        zeros(T, NR),
        zeros(T, NR),
        zeros(T, NC),
    )
end

"""
	reset!(c::Contributions) -> c

Zero every flux and source before the kernels run.
"""
function reset!(c::Contributions)
    fill!(c.F_mass, zero(eltype(c.F_mass)))
    fill!(c.F_energy, zero(eltype(c.F_energy)))
    fill!(c.S_mass, zero(eltype(c.S_mass)))
    fill!(c.S_energy, zero(eltype(c.S_energy)))
    fill!(c.θ, zero(eltype(c.θ)))
    fill!(c.r_intensity, zero(eltype(c.r_intensity)))
    fill!(c.r_velocity, zero(eltype(c.r_velocity)))
    fill!(c.F_carried, zero(eltype(c.F_carried)))
    fill!(c.S_carried, zero(eltype(c.S_carried)))
    return c
end

"""
	contribute!(c, kernel, ds, topology, material, t)

Add `kernel`'s fluxes and sources at the derived state `ds` and time `t` to `c`.
"""
function contribute! end

"""
	find_kernel(K, kernels::Tuple) -> Union{K,Nothing}

The first kernel of type `K` in `kernels`, or `nothing`; resolved at compile
time from the tuple's type.
"""
@inline find_kernel(::Type{K}, ::Tuple{}) where {K} = nothing
@inline find_kernel(::Type{K}, kernels::Tuple) where {K} =
    first(kernels) isa K ? first(kernels) : find_kernel(K, Base.tail(kernels))
