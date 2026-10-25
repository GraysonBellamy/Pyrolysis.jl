# =============================================================================
# Method of Manufactured Solutions (MMS) driver — V-MMS-01 infrastructure
# =============================================================================
#
# Closes the documented verification gap (tech ref §18.2): every other case in
# the suite compares against exact solutions whose residual source is zero.
# Here we pick smooth fields T*(z,t), ξ*ⱼ(z,t) and inject the volumetric
# sources that make them the exact solution of the FULL coupled system, which
# gives formal order-of-accuracy verification of the coupled operators.
#
# ## Continuous system targeted (Eulerian, partial-pressure transport, no radiation)
#
#   ρcp(T,ξ) ∂T/∂t = −∂z q + S_adv + Q_rxn + S_T,      q = −k_eff(T,ξ) ∂z T
#       ρcp    = Σⱼ ξⱼ cp_j(T)   (every component, gas included)
#       S_adv = −Σ_g cp_g(T) J_g ∂z T
#       Q_rxn  = −Σ_r h_r(T,ξ) r_r(ξ,T)
#   ∂ξ_g/∂t = −∂z J_g + ω̇_g + S_ξg,   J_g = −Λ(T,ξ) M_g/(RT)·∂z p_g,  p_g = ξ_g RT/(φ M_g)
#   ∂ξ_s/∂t =            ω̇_s + S_ξs                          (condensed phases)
#
# The transport reference constructs matrix weights, porosity and partial
# pressures independently from component data, then differentiates the
# continuous potential. It never calls the production transport evaluator.
#
# ## Design decisions (mirror the case docstrings before changing)
#
# 1. **Sources are injected at the state-derivative level** (K/s for T,
#    kg/(m³·s) for ξ), with the energy source divided by ρcp evaluated on
#    the *manufactured* fields — a pure function of (z,t). They enter the
#    residual through the solver's `VolumetricSources` hook, whose functions
#    of position and time the colored Jacobian differentiates with the rest
#    of the residual, so the Jacobian stays exact with MMS active.
# 2. Thermal and reaction properties use the material's declared laws. The
#    gas reference independently assembles normalized condensed matrix weights
#    and the continuous partial-pressure law, including spatially varying φ.
# 3. Boundary heat flux is positive inward; named species mass flux is positive
#    outward. Incoming enthalpy uses the manufactured boundary temperature.
#    Manufactured Neumann fluxes avoid a first-order temperature boundary
#    reconstruction entering the target truncation error.
# 4. Spatial derivatives of the manufactured fields are taken with nested
#    ForwardDiff on the field callables — no hand-derived formulas to
#    transcribe wrong. The hand-derived checks live in selftests.jl.
#
# ## ALE
# Masses are the states. Under REACTION-DRIVEN (induced) motion — the
# `MMSMotion` companion, no `w_override` — the mesh is the solid and every
# inventory rides it (the gas flux J_g is relative to it), so a z-varying
# manufactured field needs its source assembled in the material-following
# frame (`w` keyword of the source functions):
#
#   S_j = ∂tξ*_j + w·∂zξ*_j + [∂zJ*_j]_gas − ω̇*_j + ξ*_j·θ*
#
# where w·∂zξ* is the moving-point chain rule and ξ*θ* the growth of the cell
# volume at the reaction-driven strain rate θ*, with w the prefix sum of θ*
# on the manufactured fields (`mms_mesh_velocity`); the temperature source is
# assembled in the same frame. Under PRESCRIBED motion (`w_override`) the
# mesh moves over a solid at rest and every inventory, T's energy with it,
# crosses the faces at the mesh-relative velocity, so every source is the
# Eulerian one. Sources are evaluated at the cell centers.
#
# Sign conventions: q_z = −k∂zT and J positive in +z; a boundary condition's
# heat flux is positive into the domain and species mass flux positive outward.
# Thus heat uses (−q(L), +q(0)), and species use (+J(L), −J(0)).

using ForwardDiff
using LinearAlgebra: SymTridiagonal, eigen
using OrdinaryDiffEq
using Pyrolysis
using Pyrolysis.Materials: is_gas, n_components, n_reactions
using Pyrolysis: Properties, Physics

# ----------------------------------------------------------------------------
# Manufactured fields
# ----------------------------------------------------------------------------

"""
Manufactured solution fields. `T(z,t) -> Real`; `ξ` is an NTuple of
`(z,t) -> Real` callables, one per material component (same ordering as
`material.components`). All callables must be smooth (nested-ForwardDiff
differentiable) and keep ξ ≥ 0 over the case's space-time domain.
"""
struct MMSFields{FT,XT<:Tuple}
    T::FT
    ξ::XT
end

# Manufactured concentration tuple at (z,t), promoted to a homogeneous eltype
# (the property kernels require NTuple{NC,<:Real} with a single scalar type;
# a constant-field component would otherwise poison the tuple with Float64
# next to Duals).
@inline mms_xi_tuple(fields::MMSFields, z, t) = promote(map(f -> f(z, t), fields.ξ)...)

# ----------------------------------------------------------------------------
# Continuous operator evaluated on the manufactured fields
# ----------------------------------------------------------------------------

"Continuous conductive flux q(z,t) = −k_eff(T*,ξ*)·∂zT* [W/m²]."
function mms_heat_flux(material, fields::MMSFields, z, t)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    k = Properties.effective_conductivity(ξ, T, material)
    dTdz = ForwardDiff.derivative(ζ -> fields.T(ζ, t), z)
    return -k * dTdz
end

"Continuous gas mass flux J_j = −Λ M_j/(RT) ∂z p_j [kg/(m² s)]."
function mms_gas_flux(material, fields::MMSFields, j::Int, z, t)
    material.gas_transport isa LumpedTransfer ||
        error("this MMS reference requires lumped transport")
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    ids = [
        k for k in eachindex(ξ) if
        !is_gas(material.components[k]) && material.components[k].swelling_factor>0
    ]
    weights = [
        material.components[k].swelling_factor*ξ[k]/material.components[k].density(T)
        for k in ids
    ]
    weights ./= sum(weights)
    values = [material.components[k].gas_transfer_coeff(T) for k in ids]
    spec = material.mixing.λ
    parallel = sum(weights .* values)
    series =
        any(k->weights[k]>0 && values[k]==0, eachindex(weights)) ? zero(T) :
        inv(sum(weights ./ values))
    Λ =
        spec.rule==PARALLEL ? parallel :
        spec.rule==SERIES ? series :
        spec.rule==WEIGHTED ? spec.β*parallel+(1-spec.β)*series :
        error("unsupported MMS matrix closure")
    all(k->material.components[k].permeability(T)==0, ids) ||
        error("this MMS reference requires κ=0")
    potential(ζ) = begin
        Θ = fields.T(ζ, t)
        concentrations = mms_xi_tuple(fields, ζ, t)
        φ =
            1-sum(
                concentrations[k]/Pyrolysis.Materials.get_intrinsic_density(
                    material.components[k],
                    Θ,
                ) for k in eachindex(concentrations) if !is_gas(material.components[k])
            )
        φ>0 || error("manufactured fields leave the physical pore domain")
        fields.ξ[j](ζ, t)*Θ/φ
    end
    return -(Λ/T)*ForwardDiff.derivative(potential, z)
end

# Reaction rate of reaction r at manufactured state — same rate law as the
# solver (Arrhenius × ξ^n × material's own depletion limiter × tanh gates).
@inline _mms_rate(material, r::Int, ξ, T) =
    Physics.reaction_rate(material.reactions[r], ξ, T; limiter = material.depletion_limiter)

"Continuous reaction heat source Q_rxn = −Σ_r h_r·r_r [W/m³] at (z,t)."
function mms_Q_rxn(material, ξ::Tuple, T)
    NR = n_reactions(material)
    Q = zero(promote_type(eltype(ξ), typeof(T)))
    for r = 1:NR
        rate = _mms_rate(material, r, ξ, T)
        Q -= material.reactions[r].heat(T, ξ) * rate
    end
    return Q
end

"Continuous species reaction source ω̇_j [kg/(m³·s)] at (z,t)."
function mms_rxn_source(material, j::Int, ξ::Tuple, T)
    NR = n_reactions(material)
    s = zero(promote_type(eltype(ξ), typeof(T)))
    for r = 1:NR
        rxn = material.reactions[r]
        rate = _mms_rate(material, r, ξ, T)
        rxn.reactants[1][1] == j && (s -= rate)
        for (idx, ν) in rxn.products
            idx == j && (s += ν * rate)
        end
    end
    return s
end

"Continuous gas-advection energy source S_adv = −Σ_g cp_g(T)·J_g·∂zT [W/m³]."
function mms_S_adv(material, fields::MMSFields, z, t)
    T = fields.T(z, t)
    dTdz = ForwardDiff.derivative(ζ -> fields.T(ζ, t), z)
    NC = n_components(material)
    S = zero(dTdz)
    for j = 1:NC
        comp = material.components[j]
        if is_gas(comp)
            J = mms_gas_flux(material, fields, j, z, t)
            S -= comp.heat_capacity(T) * J * dTdz
        end
    end
    return S
end

"""
Manufactured energy source in **K/s** (state-derivative units — see design
decision 1):

S_T(z,t) = ∂tT* + w·∂zT* + (∂z q* − S_adv* − Q_rxn*) / ρcp(T*,ξ*)

Added to the energy row, it makes T* the exact solution of the continuous
energy equation; a function of (z,t) alone, it leaves the Jacobian exact.
`w = nothing` (fixed mesh, or prescribed motion, whose sources are Eulerian)
drops the chain-rule term; under induced motion pass the reaction-driven
`w(z,t)`.
"""
function mms_energy_source(material, fields::MMSFields, z, t; w = nothing)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    ρcp = sum(ξ[j] * material.components[j].heat_capacity(T) for j in eachindex(ξ))
    dTdt = ForwardDiff.derivative(τ -> fields.T(z, τ), t)
    div_q = ForwardDiff.derivative(ζ -> mms_heat_flux(material, fields, ζ, t), z)
    S_adv = mms_S_adv(material, fields, z, t)
    Q_rxn = mms_Q_rxn(material, ξ, T)
    S = dTdt + (div_q - S_adv - Q_rxn) / ρcp
    if w !== nothing
        S += w(z, t) * ForwardDiff.derivative(ζ -> fields.T(ζ, t), z)
    end
    return S
end

"""
Manufactured species source for component `j` in **kg/(m³·s)**:

S_ξj(z,t) = ∂tξ*_j + w·∂zξ*_j + ∂z J*_j·[gas only] − ω̇*_j + ξ*_j·θ*

`w = nothing` (fixed mesh, or prescribed motion over a solid at rest) drops
the two mesh-motion terms. Under induced motion every component rides the
mesh, so pass the reaction-driven mesh velocity `w(z,t)`
(`mms_mesh_velocity`) and the source is assembled in the material-following
frame: `w·∂zξ*` is the moving-point chain rule and `+ξ*θ*` the growth of the
cell volume at the reaction-driven strain rate θ* (`mms_theta`).
"""
function mms_species_source(material, fields::MMSFields, j::Int, z, t; w = nothing)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    dξdt = ForwardDiff.derivative(τ -> fields.ξ[j](z, τ), t)
    S = dξdt - mms_rxn_source(material, j, ξ, T)
    if is_gas(material.components[j])
        S += ForwardDiff.derivative(ζ -> mms_gas_flux(material, fields, j, ζ, t), z)
    end
    if w !== nothing
        dξdz = ForwardDiff.derivative(ζ -> fields.ξ[j](ζ, t), z)
        S += w(z, t) * dξdz + fields.ξ[j](z, t) * mms_theta(material, fields, z, t)
    end
    return S
end

# ----------------------------------------------------------------------------
# Induced (reaction-driven) mesh motion — V-MMS-04 companion machinery
# ----------------------------------------------------------------------------
#
# Under REACTION-DRIVEN motion (no `w_override`) the mesh velocity is the
# prefix sum of the volumetric strain rate θ computed from the REACTION rates
# (`Physics.volume_change_rate` on `dξ_rxn`). The manufactured counterparts
# below mirror that chain exactly — same kinetics calls, same γ/ρ weighting —
# per driver design decision 2 (never re-derive coefficients by hand).

"""
Continuous volumetric strain rate θ*(z,t) [1/s] on the manufactured fields:
θ* = Σⱼ (γⱼ/ρⱼ(T*))·dξⱼ_rxn — REACTION rates only, via the same
`Physics.volume_change_rate` the mesh closure's strain rates use.
"""
function mms_theta(material, fields::MMSFields, z, t)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    NC = n_components(material)
    dξ_rxn = ntuple(j -> mms_rxn_source(material, j, ξ, T), NC)
    return Physics.volume_change_rate(promote(dξ_rxn...), T, material)
end

# Golub–Welsch Gauss–Legendre nodes/weights on [−1, 1].
function _gauss_legendre(n::Int)
    β = [k / sqrt(4.0 * k^2 - 1.0) for k = 1:(n-1)]
    x, V = eigen(SymTridiagonal(zeros(n), β))
    return x, 2 .* abs2.(V[1, :])
end

const _MMS_GL_NODES, _MMS_GL_WEIGHTS = _gauss_legendre(32)

"""
Continuous mesh velocity under reaction-driven motion:
w*(z,t) = ∫₀ᶻ θ*(ζ,t) dζ, by 32-node Gauss–Legendre quadrature — machine
precision for the smooth θ* fields used in MMS cases, with the package
kinetics as the θ source of truth.
"""
function mms_mesh_velocity(material, fields::MMSFields, z, t)
    half = z / 2
    s = zero(half)
    @inbounds for k in eachindex(_MMS_GL_NODES)
        ζ = half * (_MMS_GL_NODES[k] + 1)
        s += _MMS_GL_WEIGHTS[k] * mms_theta(material, fields, ζ, t)
    end
    return half * s
end

"""
Opt-in companion for MMS cases with REACTION-DRIVEN (induced) mesh motion.
`w(z,t)` must equal ∫₀ᶻ θ* dζ on the manufactured fields — build it with
[`mms_mesh_velocity`](@ref) (or a case-local closed form verified against
it). Passed to `solve_mms(...; motion=...)`, it switches every field's
manufactured source to the material-following frame (`w` keyword of
`mms_energy_source` / `mms_species_source`).
"""
struct MMSMotion{WF}
    w::WF
end

# ----------------------------------------------------------------------------
# Manufactured boundary conditions
# ----------------------------------------------------------------------------

"""
Time-dependent `HeatFluxBC` pair consistent with the manufactured field on
z ∈ [0, L]. Into-domain convention: top = −q*(L,t), bottom = +q*(0,t).
"""
function mms_thermal_bcs(material, fields::MMSFields, L)
    return (
        top = HeatFluxBC(t -> -mms_heat_flux(material, fields, L, t)),
        bottom = HeatFluxBC(t -> mms_heat_flux(material, fields, 0.0, t)),
    )
end

"Named outward species fluxes with complete incoming thermodynamic data."
function mms_mass_bcs(material, fields::MMSFields, L)
    ids = findall(is_gas, collect(material.components))
    names = Tuple(material.components[j].name for j in ids)
    boundary(z, sign) = MassFluxBC(
        fluxes = NamedTuple{names}(
            Tuple(t->sign*mms_gas_flux(material, fields, j, z, t) for j in ids),
        ),
        inflow_temperature = NamedTuple{names}(Tuple(t->fields.T(z, t) for j in ids)),
    )
    return (top = boundary(L, 1), bottom = boundary(0.0, -1))
end

# ----------------------------------------------------------------------------
# Solve path — the solver's volumetric-source hook
# ----------------------------------------------------------------------------

"""
Minimal solution container for MMS runs (the harness only needs `retcode`
plus state trajectories).

Fields: `t` (n_times), `T` (n_cells × n_times), `ξ` (NC × n_cells × n_times),
`z` (n_cells × n_times cell centers — constant columns on a fixed mesh), and
`statistics`, the solve's `SolverStatistics`.
"""
struct MMSSolution
    retcode::Symbol
    t::Vector{Float64}
    T::Matrix{Float64}
    ξ::Array{Float64,3}
    z::Matrix{Float64}
    statistics::Any
end

"""
    mms_top_node(sol, j) -> z_top

Top node position at snapshot `j`, reconstructed from the saved cell centers:
the bottom node is pinned at 0 and cᵢ = (zᵢ + zᵢ₊₁)/2, so the node ladder
zᵢ₊₁ = 2cᵢ − zᵢ is exact.
"""
function mms_top_node(sol::MMSSolution, j::Int)
    z = 0.0
    for i in axes(sol.z, 1)
        z = 2.0 * sol.z[i, j] - z
    end
    return z
end

"""
    solve_mms(problem, fields; kwargs...) -> MMSSolution

Solve `problem` with the manufactured sources for `fields`, entering through
the solver's `VolumetricSources` hook in the energy state `energy_state`: the
species sources on the mass rows as `S_ξ V`, the temperature-rate source on
the energy row (multiplied by the cell's heat capacity in the inventory
form). Both are functions of `(z, t)` at the cell centers, which the colored
Jacobian differentiates with the rest of the residual.

Keywords:
- `use_ale=false`: move the mesh.
- `w_override=nothing`: prescribed mesh-velocity function `(z_node, t) -> w`
  (ALE only), moving the mesh over a solid at rest; every source is then the
  Eulerian one.
- `motion=nothing`: opt-in [`MMSMotion`](@ref) for REACTION-DRIVEN mesh
  motion (ALE only, no `w_override`): every source — T and all species —
  switches to the material-following frame with `motion.w`.
- `dt=nothing`: fixed step (sets `adaptive=false`) for temporal-order sweeps.
- `saveat=nothing`: defaults to `problem.output_times` if non-empty, else the
  final time.
"""
function solve_mms(
    problem,
    fields::MMSFields;
    use_ale::Bool = false,
    w_override = nothing,
    motion::Union{Nothing,MMSMotion} = nothing,
    abstol::Real = 1e-10,
    reltol::Real = 1e-8,
    dt::Union{Nothing,Real} = nothing,
    saveat = nothing,
    maxiters::Int = 10^6,
    energy_state::Symbol = :temperature,
)
    material = problem.material
    save_times = if saveat !== nothing
        saveat
    elseif !isempty(problem.output_times)
        problem.output_times
    else
        [problem.tspan[2]]
    end
    (w_override !== nothing || motion !== nothing) &&
        !use_ale &&
        error("w_override and motion require use_ale = true")
    w_override !== nothing &&
        motion !== nothing &&
        error("w_override (prescribed motion) and motion (induced) are exclusive")
    NC = n_components(material)
    w = motion === nothing ? nothing : motion.w
    mass = (z, t) -> ntuple(j -> mms_species_source(material, fields, j, z, t; w), Val(NC))
    energy = (z, t) -> mms_energy_source(material, fields, z, t; w)
    step_kwargs = dt === nothing ? (;) : (dt = dt, adaptive = false)
    sol = solve(
        problem;
        energy_state,
        use_ale,
        w_override,
        abstol,
        reltol,
        saveat = save_times,
        maxiters,
        sources = Pyrolysis.Assembly.VolumetricSources(mass, energy),
        step_kwargs...,
    )
    return MMSSolution(sol.retcode, sol.t, sol.T, sol.ξ, sol.cell_centers, sol.statistics)
end

# ----------------------------------------------------------------------------
# Error norms against the manufactured fields
# ----------------------------------------------------------------------------

"L∞ temperature error at snapshot index `m` (default: final time)."
function mms_Linf_T(sol::MMSSolution, fields::MMSFields; m::Int = length(sol.t))
    t = sol.t[m]
    return maximum(abs(sol.T[i, m] - fields.T(sol.z[i, m], t)) for i in axes(sol.T, 1))
end

"Cell-averaged L₂ temperature error at snapshot index `m` (default: final)."
function mms_L2_T(sol::MMSSolution, fields::MMSFields; m::Int = length(sol.t))
    t = sol.t[m]
    n = size(sol.T, 1)
    return sqrt(sum(abs2(sol.T[i, m] - fields.T(sol.z[i, m], t)) for i = 1:n) / n)
end

"L∞ error of component `j`'s concentration at snapshot `m` (default: final)."
function mms_Linf_xi(sol::MMSSolution, fields::MMSFields, j::Int; m::Int = length(sol.t))
    t = sol.t[m]
    return maximum(
        abs(sol.ξ[j, i, m] - fields.ξ[j](sol.z[i, m], t)) for i in axes(sol.T, 1)
    )
end
