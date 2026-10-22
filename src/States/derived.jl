# ============================================================================
# Per-cell reference data
# ============================================================================

"""
	CellReference

Constants fixed at setup, per cell unless stated:

- `volume`: initial volume `Vᵢ(0)` [m³]
- `solid_mass`: the scale `mˢᵢ(0)` [kg] of the cell's mass tolerances: its
  initial condensed (solid and liquid) mass, or for a cell that starts with
  none its total initial mass, or 1 kg/m³ times its volume
- `mass_scale`: the positive-part scale `δᵢ` [kg] of the caloric extension and
  of the property evaluation
- `volume_floor`: smooth floor [m³] on the volume used for `m/V`
- `capacity_floor`: smooth floor [J/K] on the cell's heat-capacity sum
- `positions`: reference node positions [m], the geometry of a fixed mesh,
  the `D` coordinates of every node in turn
- `thickness`: initial size `hᵢ(0)` [m]: on a line the distance between the
  cell's two nodes, in more dimensions `Vᵢ(0)^(1/D)`
- `cross_section`: the cross-section measure `A₀` of the mesh
"""
struct CellReference
    volume::Vector{Float64}
    solid_mass::Vector{Float64}
    mass_scale::Vector{Float64}
    volume_floor::Vector{Float64}
    capacity_floor::Vector{Float64}
    positions::Vector{Float64}
    thickness::Vector{Float64}
    cross_section::Float64
end

# ============================================================================
# Gas transport model
# ============================================================================

"""
	GasTransport

Which pore-gas properties the state mapping evaluates for gas transport.
"""
abstract type GasTransport end

"""
	NoGasTransport <: GasTransport

No gas moves: no pore-gas property is evaluated.
"""
struct NoGasTransport <: GasTransport end

"""
	FickTransport <: GasTransport

Diffusive transport `J = −(λ/T) ∂(ξT)/∂n`: the cell gas-transfer coefficients.
"""
struct FickTransport <: GasTransport end

"""
	DarcyFickTransport <: GasTransport

Diffusive and pressure-driven transport: the gas-transfer coefficients, the
porosity, the ideal-gas pressure and the Darcy mobility `κ/μ`.
"""
struct DarcyFickTransport <: GasTransport end

# Smooth floor on the porosity the pressure divides by, the Darcy closure's
# floor on the gas volume fraction.
const POROSITY_FLOOR = 1e-12

# Smooth floor on the cell conductivity [W/(m·K)], far below any material's:
# a trial state with no mass left has none, and the face temperatures of the
# boundary conditions divide by it.
const CONDUCTIVITY_FLOOR = 1e-6

# ============================================================================
# State setup
# ============================================================================

"""
	StateSetup

Everything the state mapping reads besides the state vector: the mesh
topology, the material, its caloric laws, the per-cell reference data, for
the inventory form the caloric extension (`nothing` for the temperature
form), and the [`GasTransport`](@ref) model.
"""
struct StateSetup{TopT<:MeshTopology,Mat<:Material,LawT<:CaloricLaws,CalT,G<:GasTransport}
    topology::TopT
    material::Mat
    laws::LawT
    reference::CellReference
    caloric::CalT
    gas::G
end

"""
	StateSetup(topology, material, reference, caloric; gas = NoGasTransport())

The setup with the material's caloric laws built.
"""
StateSetup(
    topo::MeshTopology,
    material::Material,
    reference::CellReference,
    caloric;
    gas::GasTransport = NoGasTransport(),
) = StateSetup(topo, material, CaloricLaws(material), reference, caloric, gas)

# ============================================================================
# Derived state
# ============================================================================

"""
	DerivedState{NC,T}

The quantities every kernel reads, mapped from the state vector by
[`derive!`](@ref) with the state's element type `T` (a dual number under
automatic differentiation):

- `m`, `ξ`, `ξ⁺`: component masses (as stored), concentrations `m/V`, and the
  positive-part concentrations the property laws evaluate at (`NC × n`)
- `T`, `T_trial`: cell temperatures and their smooth clamp into the property
  bracket
- `V`, `V_eval`: signed geometric cell volumes and their smooth floor
- `cp_sum`: the smoothly floored heat-capacity sum `Σⱼ mⱼ cⱼ(T_trial)` [J/K]
- `h`: component sensible enthalpies at the cell temperature (`NC × n`)
- `x`, `w`, `A`: node positions, node velocities (moving mesh) and the
  cross-section measure; positions and velocities hold the `D` components of
  every node in turn
- `multiplier`: the closure's multipliers per cell (moving mesh,
  multipliers), from which the closure fills the node velocities
- `c`: the node corrections `c = w − vₛ`, the mesh velocity relative to the
  solid's (zero without a correction), filled by `Closures.correction!`
- `θ_A`: the lateral strain rate of the cross section [1/s] (variable cross
  section; zero otherwise)
- `G`, `M`: the conversion-history pair per cell [kg] (history pair)
- `F`: the deformation gradient per cell, column by column (`D² × n`,
  deformation gradient)
- `geometry`: the mesh metrics at `x`
- `k`: cell thermal conductivities, smoothly floored at `CONDUCTIVITY_FLOOR`
- `ε`, `a`: cell surface emissivities and absorptivities (the mixing-weighted
  sums, which stay in `[0, 1]`)
- `I`, `α`: face intensities [W/m²] as stored, and cell absorption
  coefficients at the positive-part concentrations [1/m] (in-depth radiation)
- `λ`: cell gas-transfer coefficients [m²/s] (gas transport)
- `φ`, `P`, `mob`: the smoothly floored porosity, the ideal-gas pressure [Pa]
  and the Darcy mobility `κ/μ` with a smoothly floored viscosity [m²/(Pa·s)]
  (Darcy–Fick transport)
- `active`: the live cells, those whose rows the layout's mask keeps; the
  kernels skip every face and patch beyond them
- `T_warm`: the latest temperatures of the plain (non-dual) mapping, the
  starting point of the inventory inversion
"""
struct DerivedState{NC,T,G<:MeshGeometry}
    m::Matrix{T}
    ξ::Matrix{T}
    ξ⁺::Matrix{T}
    T::Vector{T}
    T_trial::Vector{T}
    V::Vector{T}
    V_eval::Vector{T}
    cp_sum::Vector{T}
    h::Matrix{T}
    x::Vector{T}
    w::Vector{T}
    multiplier::Vector{T}
    c::Vector{T}
    A::Base.RefValue{T}
    θ_A::Base.RefValue{T}
    G::Vector{T}
    M::Vector{T}
    F::Matrix{T}
    geometry::G
    k::Vector{T}
    ε::Vector{T}
    a::Vector{T}
    I::Vector{T}
    α::Vector{T}
    λ::Vector{T}
    φ::Vector{T}
    P::Vector{T}
    mob::Vector{T}
    active::BitVector
    T_warm::Vector{Float64}
end

"""
	DerivedState{T}(layout, setup, T_warm) -> DerivedState

Buffers with element type `T` for `layout`, the geometry of a fixed mesh
already filled. `T_warm` is shared by every derived state of one problem.
"""
function DerivedState{T}(
    layout::StateLayout{NC},
    setup::StateSetup,
    T_warm::Vector{Float64},
) where {T,NC}
    topo = setup.topology
    n = n_cells(layout)
    nn = dimension(layout) * n_nodes(layout)
    ref = setup.reference
    geometry = MeshGeometry{T}(topo)
    x = T.(ref.positions)
    A = Ref{T}(T(ref.cross_section))
    geometry!(geometry, x, A[], topo)
    return DerivedState{NC,T,typeof(geometry)}(
        zeros(T, NC, n),
        zeros(T, NC, n),
        zeros(T, NC, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, NC, n),
        x,
        zeros(T, nn),
        zeros(T, n),
        zeros(T, nn),
        A,
        Ref{T}(zero(T)),
        zeros(T, n),
        zeros(T, n),
        zeros(T, dimension(layout)^2, n),
        geometry,
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n_faces(topo)),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n),
        trues(n),
        T_warm,
    )
end

"""
	derive!(ds, u, layout, setup) -> ds

Map the state vector `u` to the derived state: the only reader of `u`.

Order: the live cells from the layout's mask; the geometry from the node
positions (a fixed mesh keeps its precomputed geometry), with the node
velocities or the closure's multipliers, and under a variable cross section
the area and its lateral rate; the carried fields; the face intensities; the
masses and their concentrations, the cell temperatures (the stored
temperature, or the inversion of the enthalpy inventory), their trial clamp,
the heat-capacity sum, the sensible enthalpies, the cell properties with the
absorption coefficients under in-depth radiation, and the pore-gas properties
the gas transport model needs.
"""
function derive!(
    ds::DerivedState{NC,T},
    u::AbstractVector,
    layout::StateLayout{NC},
    setup::StateSetup,
) where {NC,T}
    topo = setup.topology
    ref = setup.reference
    material = setup.material
    n = n_cells(layout)

    @inbounds for i = 1:n
        ds.active[i] = layout.mask[energy_index(layout, i)]
    end
    if is_moving(layout)
        _motion!(ds, u, layout)
        if has_chi(layout)
            ds.A[] = u[area_index(layout)]
            ds.θ_A[] = u[lateral_rate_index(layout)]
        end
        geometry!(ds.geometry, ds.x, ds.A[], topo)
    end
    if has_history(layout)
        @inbounds for i = 1:n
            ds.G[i] = u[history_index(layout, i, 1)]
            ds.M[i] = u[history_index(layout, i, 2)]
        end
    end
    if has_deformation(layout)
        D = dimension(layout)
        @inbounds for i = 1:n, b = 1:D, a = 1:D
            ds.F[D*(b-1)+a, i] = u[deformation_index(layout, i, a, b)]
        end
    end
    if has_intensities(layout)
        @inbounds for f in eachindex(ds.I)
            ds.I[f] = u[intensity_index(layout, f)]
        end
    end

    @inbounds for i = 1:n
        V = ds.geometry.cell_volume[i]
        Vf = smooth_floor(V, ref.volume_floor[i])
        ds.V[i] = V
        ds.V_eval[i] = Vf
        δ = ref.mass_scale[i]
        for j = 1:NC
            m = u[mass_index(layout, i, j)]
            ds.m[j, i] = m
            ds.ξ[j, i] = m / Vf
            ds.ξ⁺[j, i] = positive_part(m, δ) / Vf
        end
    end

    _temperatures!(ds, u, layout, setup, energy_state(layout))

    @inbounds for i = 1:n
        Tt = trial_temperature(ds.T[i])
        ds.T_trial[i] = Tt
        m = ntuple(j -> ds.m[j, i], Val(NC))
        C = heat_capacity_sum(zero(T), setup.laws, m, Tt)
        ds.cp_sum[i] = smooth_floor(C, ref.capacity_floor[i])
        h = component_enthalpies(setup.laws, ds.T[i])
        for j = 1:NC
            ds.h[j, i] = h[j]
        end
        ξ⁺ = ntuple(j -> ds.ξ⁺[j, i], Val(NC))
        k = effective_conductivity(ξ⁺, Tt, material; skip_absent = Val(true))
        ds.k[i] = smooth_floor(k, CONDUCTIVITY_FLOOR)
        ds.ε[i] = mixture_emissivity(ξ⁺, Tt, material)
        ds.a[i] = mixture_absorptivity(ξ⁺, Tt, material)
        _gas_properties!(ds, i, ξ⁺, Tt, material, setup.gas)
        has_intensities(layout) && (ds.α[i] = effective_absorption(ξ⁺, material))
    end
    return ds
end

# The node positions and the closure's unknowns of a moving mesh.
@inline function _motion!(ds, u, layout::StateLayout)
    D = dimension(layout)
    @inbounds for nd = 1:n_nodes(layout), k = 1:D
        ds.x[D*(nd-1)+k] = u[position_index(layout, nd, k)]
    end
    if closure_unknowns(layout) === :velocities
        @inbounds for nd = 1:n_nodes(layout), k = 1:D
            ds.w[D*(nd-1)+k] = u[velocity_index(layout, nd, k)]
        end
    else
        @inbounds for i = 1:n_cells(layout)
            ds.multiplier[i] = u[multiplier_index(layout, i)]
        end
    end
    return nothing
end

# Pore-gas properties of cell `i` at the positive-part concentrations and the
# trial temperature.
@inline _gas_properties!(ds, i, ξ⁺, Tt, material, ::NoGasTransport) = nothing

@inline function _gas_properties!(ds, i, ξ⁺, Tt, material, ::FickTransport)
    ds.λ[i] = effective_gas_transfer(ξ⁺, Tt, material)
    return nothing
end

@inline function _gas_properties!(ds, i, ξ⁺, Tt, material, ::DarcyFickTransport)
    ds.λ[i] = effective_gas_transfer(ξ⁺, Tt, material)
    φ = smooth_floor(1 - skeletal_volume_fraction(ξ⁺, Tt, material), POROSITY_FLOOR)
    ds.φ[i] = φ
    ds.P[i] = ideal_gas_pressure(ξ⁺, Tt, φ, material)
    μ = smooth_floor(material.gas_viscosity(Tt), MU_FLOOR)
    ds.mob[i] = effective_permeability(ξ⁺, Tt, material, φ) / μ
    return nothing
end

@inline function _temperatures!(
    ds::DerivedState,
    u::AbstractVector,
    layout::StateLayout,
    setup::StateSetup,
    ::TemperatureState,
)
    @inbounds for i = 1:n_cells(layout)
        ds.T[i] = u[energy_index(layout, i)]
    end
    return nothing
end

function _temperatures!(
    ds::DerivedState{NC,T},
    u::AbstractVector,
    layout::StateLayout,
    setup::StateSetup,
    ::EnthalpyState,
) where {NC,T}
    ext = setup.caloric
    laws = setup.laws
    @inbounds for i = 1:n_cells(layout)
        m = ntuple(j -> ds.m[j, i], Val(NC))
        Ti, _ = invert_caloric(ext, laws, u[energy_index(layout, i)], m, i, ds.T_warm[i])
        ds.T[i] = Ti
        T === Float64 && isfinite(Ti) && (ds.T_warm[i] = Ti)
    end
    return nothing
end

"""
	energy_state_value(energy, setup, T, m, i)

The energy state of cell `i` with temperature `T` and masses `m`: `T` itself in
the temperature form; in the inventory form the extended relation of the
cell's [`CaloricExtension`](@ref), so the inversion returns `T` exactly.
"""
@inline energy_state_value(::TemperatureState, setup, T, m, i) = T
@inline energy_state_value(::EnthalpyState, setup, T, m, i) =
    extended_energy(setup.caloric, setup.laws, T, m, i)
