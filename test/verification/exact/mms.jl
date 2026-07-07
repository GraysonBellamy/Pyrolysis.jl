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
# ## Continuous system targeted (Eulerian, Fickian transport, no radiation,
# ## energy_form = :standard — the discretization's consistency target)
#
#   ρcp_eff(T,ξ) ∂T/∂t = −∂z q + S_conv + Q_rxn + S_T,      q = −k_eff(T,ξ) ∂z T
#       S_conv = −Σ_g cp_g(T) J_g ∂z T
#       Q_rxn  = −Σ_r h_r(T,ξ) r_r(ξ,T)
#   ∂ξ_g/∂t = −∂z J_g + ω̇_g + S_ξg,   J_g = −(λ_eff(T,ξ)/T)·∂z(ξ_g T)
#   ∂ξ_s/∂t =            ω̇_s + S_ξs                          (condensed phases)
#
# The gas flux is the ThermaKin volume-fraction-gradient form (see
# src/Physics/gas_transport.jl), NOT plain Fick — the manufactured source must
# match it exactly or the "error" never converges.
#
# ## Design decisions (mirror the case docstrings before changing)
#
# 1. **Sources are injected at the state-derivative level** (K/s for T,
#    kg/(m³·s) for ξ), with the energy source divided by ρcp_eff evaluated on
#    the *manufactured* fields — a pure function of (z,t). A W/m³ source
#    divided by the numerical ρcp(u) would be state-dependent and silently
#    invalidate the analytic Jacobian; this way ∂S/∂u ≡ 0 and the shipped
#    `Structured` Jacobian remains exact with MMS active.
# 2. **Every coefficient is evaluated through the package's own property and
#    kinetics functions** (`effective_heat_capacity`, `effective_conductivity`,
#    `effective_gas_transfer`, `reaction_rate` with the material's own
#    depletion limiter). The continuous operator therefore matches the
#    discretization's consistency target by construction — including mixing
#    rules, floors, tanh gates and limiter — with no re-derivation to drift.
# 3. **Boundary conditions are manufactured Neumann fluxes**
#    (`HeatFluxBC(t->…)` / `MassFluxBC(t->…)`, positive into domain): the FV
#    scheme imposes face fluxes exactly, so the boundary contributes no
#    order-reducing reconstruction (the Dirichlet path is first-order at the
#    face, see finite_volume.jl `dz_half` reconstruction).
# 4. Spatial derivatives of the manufactured fields are taken with nested
#    ForwardDiff on the field callables — no hand-derived formulas to
#    transcribe wrong. The hand-derived checks live in selftests.jl.
#
# ## ALE
# The same injection works on a moving mesh for T and GAS fields: on a grid
# point z_i(t) moving with mesh velocity w, d/dt T*(z_i(t),t) = ∂tT* + w ∂zT*,
# and the discrete ALE advection term supplies (its approximation of) w ∂zT*,
# so the Eulerian source formula evaluated at the CURRENT cell center is
# exactly right — the wrapper reads cell centers from the live z-block. Mesh
# motion is prescribed through the `Workspace.w_override` hook (free-stream &
# ALE-MMS cases).
#
# CONSTRAINT — condensed phases are mesh-bound: the ALE formulation tracks
# solids/liquids Lagrangian-ly (v_solid = w ⇒ NO advection term; see
# compute_ale_species_advection!), so under a PRESCRIBED motion that is not
# the material velocity, a z-varying manufactured solid field has no
# discrete operator to balance its w·∂zξ* term and the error cannot
# converge. ALE-MMS cases must keep every condensed-phase manufactured field
# spatially uniform (time variation is fine — it still exercises the source
# path); gas fields and T carry the full space-time structure.
#
# Under REACTION-DRIVEN (induced) motion — no w_override — the mesh IS the
# solid material frame, so z-varying condensed fields become representable;
# their manufactured sources must then be assembled in the material-following
# frame (S_s = ∂tξ* + w*·∂zξ* − ω̇* + ξ*·θ*). That path is the opt-in
# `MMSMotion` companion below (V-MMS-04).
#
# Sign conventions (verified against finite_volume.jl):
#   q_z = −k∂zT, J positive in +z; "into domain" flux at top = −(flux at L),
#   at bottom = +(flux at 0). Mass BCs: J[j,f] += σ·flux_into_domain with
#   σ = −1 (top), +1 (bottom); thermal BCs write q[f] = ∓q_conducted.

using ForwardDiff
using LinearAlgebra: SymTridiagonal, eigen
using OrdinaryDiffEq
using Pyrolysis
using Pyrolysis.Materials: is_gas, n_components, n_reactions
using Pyrolysis: Properties, Physics
# The default integrator lives in the SDIRK sub-library, not the umbrella
# (same split solve() works around — see src/Solver/Solver.jl). AutoFiniteDiff
# is needed for the :fd Jacobian path: the canonical workspace is Float64-typed,
# so the integrator's default ForwardDiff cannot pass Duals through the wrapped
# RHS (solve() makes the same choice for its no-Jacobian fallback).
using Pyrolysis.Solver: KenCarp4, AutoFiniteDiff

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

"""
Continuous gas mass flux of component `j` [kg/(m²·s)]:
J_j = −(λ_eff/T)·∂z(ξ_j T) — the ThermaKin volume-fraction-gradient form the
discrete `gas_flux` kernel approximates.
"""
function mms_gas_flux(material, fields::MMSFields, j::Int, z, t)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    λ = Properties.effective_gas_transfer(ξ, T, material)
    dξT = ForwardDiff.derivative(ζ -> fields.ξ[j](ζ, t) * fields.T(ζ, t), z)
    return -(λ / T) * dξT
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

"Continuous gas-advection energy source S_conv = −Σ_g cp_g(T)·J_g·∂zT [W/m³]."
function mms_S_conv(material, fields::MMSFields, z, t)
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

S_T(z,t) = ∂tT* + (∂z q* − S_conv* − Q_rxn*) / ρcp_eff(T*,ξ*)

Added to `du` on the temperature block, it makes T* the exact solution of the
continuous energy equation while leaving the analytic Jacobian untouched.
"""
function mms_energy_source(material, fields::MMSFields, z, t)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    ρcp = Properties.effective_heat_capacity(ξ, T, material)
    dTdt = ForwardDiff.derivative(τ -> fields.T(z, τ), t)
    div_q = ForwardDiff.derivative(ζ -> mms_heat_flux(material, fields, ζ, t), z)
    S_conv = mms_S_conv(material, fields, z, t)
    Q_rxn = mms_Q_rxn(material, ξ, T)
    return dTdt + (div_q - S_conv - Q_rxn) / ρcp
end

"""
Manufactured species source for component `j` in **kg/(m³·s)**:

S_ξj(z,t) = ∂tξ*_j + ∂z J*_j·[gas only] − ω̇*_j
"""
function mms_species_source(material, fields::MMSFields, j::Int, z, t)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    dξdt = ForwardDiff.derivative(τ -> fields.ξ[j](z, τ), t)
    src = mms_rxn_source(material, j, ξ, T)
    if is_gas(material.components[j])
        divJ = ForwardDiff.derivative(ζ -> mms_gas_flux(material, fields, j, ζ, t), z)
        return dξdt + divJ - src
    else
        return dξdt - src
    end
end

# ----------------------------------------------------------------------------
# Induced (reaction-driven) mesh motion — V-MMS-04 companion machinery
# ----------------------------------------------------------------------------
#
# Under REACTION-DRIVEN motion (no `w_override`) the mesh velocity is the
# prefix sum of the volumetric strain rate θ computed from the REACTION rates
# (`_compute_ale_mesh_velocity!` → `Physics.volume_change_rate` on
# `dξ_rxn`). The manufactured counterparts below mirror that chain exactly —
# same kinetics calls, same γ/ρ weighting — per driver design decision 2
# (never re-derive coefficients by hand).
#
# Condensed phases are mesh-bound (no advection operator; dilation −ξθ
# instead), so with the mesh moving at the SOLID material velocity the
# manufactured condensed-phase source must be assembled in the
# material-following frame:
#
#   S_s = ∂tξ* + w*·∂zξ* − ω̇* + ξ*·θ*   (evaluated at live cell centers)
#
# — the `w*·∂zξ*` term is the moving-point chain rule (no discrete operator
# supplies it for solids) and `+ξ*·θ*` balances the discrete dilation
# correction. Gas and T sources keep the existing Eulerian forms. All of this
# is OPT-IN through the `MMSMotion` companion so V-MMS-01/02/03 are untouched.

"""
Continuous volumetric strain rate θ*(z,t) [1/s] on the manufactured fields:
θ* = Σⱼ (γⱼ/ρⱼ(T*))·dξⱼ_rxn — REACTION rates only, via the same
`Physics.volume_change_rate` the residual's ALE mesh velocity uses.
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
precision for the smooth θ* fields used in MMS cases (cross-cutting
requirement 4), with the package kinetics as the θ source of truth.
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
it). Passed to `solve_mms(...; motion=...)`, it switches every condensed
(non-gas) component's manufactured source to the material-following frame
(`mms_species_source_material`); gas and T sources keep the Eulerian forms.
"""
struct MMSMotion{WF}
    w::WF
end

"""
Manufactured condensed-phase source in the material-following frame
[kg/(m³·s)] (see the section comment above):

S_j(z,t) = ∂tξ*_j + w*·∂zξ*_j − ω̇*_j + ξ*_j·θ*
"""
function mms_species_source_material(
    material,
    fields::MMSFields,
    motion::MMSMotion,
    j::Int,
    z,
    t,
)
    T = fields.T(z, t)
    ξ = mms_xi_tuple(fields, z, t)
    dξdt = ForwardDiff.derivative(τ -> fields.ξ[j](z, τ), t)
    dξdz = ForwardDiff.derivative(ζ -> fields.ξ[j](ζ, t), z)
    ω̇ = mms_rxn_source(material, j, ξ, T)
    θ = mms_theta(material, fields, z, t)
    return dξdt + motion.w(z, t) * dξdz - ω̇ + fields.ξ[j](z, t) * θ
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

"""
Time-dependent `MassFluxBC` pair for the manufactured gas flux. `MassFluxBC`
applies the same flux to EVERY gas component (see `mass_bc_flux`), so the
material must have exactly one gas component — asserted here.
"""
function mms_mass_bcs(material, fields::MMSFields, L)
    gas_idx = findall(is_gas, collect(material.components))
    length(gas_idx) == 1 || error(
        "mms_mass_bcs requires exactly one gas component (MassFluxBC applies " *
        "one flux to every gas species); got $(length(gas_idx))",
    )
    j = gas_idx[1]
    return (
        top = MassFluxBC(t -> -mms_gas_flux(material, fields, j, L, t)),
        bottom = MassFluxBC(t -> mms_gas_flux(material, fields, j, 0.0, t)),
    )
end

# ----------------------------------------------------------------------------
# Solve path — wrapped-RHS ODEProblem around the pure residual layer
# ----------------------------------------------------------------------------

"""
Minimal solution container for MMS runs (the harness only needs `retcode`
plus state trajectories; `PyrolysisSolution`'s extras don't apply to the
wrapped-RHS path).

Fields: `t` (n_times), `T` (n_cells × n_times), `ξ` (NC × n_cells × n_times),
`z` (n_cells × n_times cell centers — constant columns for Eulerian).
"""
struct MMSSolution
    retcode::Symbol
    t::Vector{Float64}
    T::Matrix{Float64}
    ξ::Array{Float64,3}
    z::Matrix{Float64}
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

Solve `problem` with the manufactured sources for `fields` injected into the
pure residual (`residual!`). Replicates `solve()`'s ODE setup (same default
integrator and analytic-Jacobian attachment) but wraps the RHS — `solve()`
rejects unknown kwargs and has no volumetric source hook, so this IS the
supported injection point (design note in the file header).

Keywords:
- `use_ale=false`: ALE mode (cell centers for source evaluation are then read
  from the live z-block each RHS call).
- `w_override=nothing`: prescribed mesh-velocity function `(z_node, t) -> w`,
  forwarded to the workspace hook (ALE only; requires `jacobian = :fd` — the
  structured geometry Jacobian assumes reaction-driven motion).
- `motion=nothing`: opt-in [`MMSMotion`](@ref) for REACTION-DRIVEN mesh
  motion (ALE only, no `w_override`): condensed-phase sources switch to the
  material-following frame (`mms_species_source_material`); gas and T keep
  the Eulerian forms.
- `jacobian=:analytic`: `:analytic` attaches the `Structured(DirectSolve())`
  Jacobian; `:approx` attaches `Structured(ApproxSolve())` (valid only for
  reaction-driven motion); `:fd` lets the integrator use its own colored
  finite differences. NOTE (induced-motion caveat): with sources evaluated at
  live cell centers, `du += S(z(u), t)` has a real ∂S/∂z chain that neither
  structured backend carries — harmless for solution correctness and observed
  order (inexact Newton matrix), and DirectSolve-vs-ApproxSolve comparisons
  stay clean because both share the omission.
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
    jacobian::Symbol = :analytic,
    integrator = KenCarp4(),
    abstol::Real = 1e-10,
    reltol::Real = 1e-8,
    dt::Union{Nothing,Real} = nothing,
    saveat = nothing,
    maxiters::Int = 10^6,
)
    material = problem.material
    # The problem's callable ICs live on the PyrolysisProblem, not the mesh;
    # initialize_state applies them before packing u0 (same path as solve()).
    def, ws, u0 = initialize_state(problem; use_ale = use_ale)
    if w_override !== nothing
        use_ale || error("w_override requires use_ale = true")
        motion === nothing ||
            error("w_override (prescribed motion) and motion (induced) are exclusive")
        jacobian === :fd || error(
            "w_override requires jacobian = :fd: the structured geometry " *
            "Jacobian assumes reaction-driven mesh velocity (see Workspace.w_override)",
        )
        ws.w_override = w_override
    end
    motion === nothing ||
        use_ale ||
        error("motion (induced mesh motion) requires use_ale = true")

    n = ws.layout.n_cells
    NC = n_components(material)
    zoff = n * (1 + NC)
    zc0 = copy(ws.state_cache.cell_centers)
    gas_mask = ntuple(j -> is_gas(material.components[j]), NC)

    rhs! = function (du, u, p, t)
        residual!(du, u, def, t, ws)
        @inbounds for i = 1:n
            z_i = use_ale ? (u[zoff+i] + u[zoff+i+1]) / 2 : zc0[i]
            du[i] += mms_energy_source(material, fields, z_i, t)
            for j = 1:NC
                du[n+(j-1)*n+i] += if motion !== nothing && !gas_mask[j]
                    mms_species_source_material(material, fields, motion, j, z_i, t)
                else
                    mms_species_source(material, fields, j, z_i, t)
                end
            end
        end
        return nothing
    end

    tspan = problem.tspan
    ode_problem = if jacobian === :analytic || jacobian === :approx
        strategy = jacobian === :approx ? ApproxSolve() : DirectSolve()
        spec = JacobianSpec(Structured(strategy = strategy), def, ws)
        jac! = (J, u, p, t) -> compute_jacobian!(J, spec, u, def, t, ws)
        proto = if spec.sparsity_phys !== nothing
            S =
                spec.sparsity_geom === nothing ? spec.sparsity_phys :
                (spec.sparsity_phys .| spec.sparsity_geom)
            Float64.(S)
        else
            nothing
        end
        fn =
            proto === nothing ? OrdinaryDiffEq.ODEFunction(rhs!; jac = jac!) :
            OrdinaryDiffEq.ODEFunction(rhs!; jac = jac!, jac_prototype = proto)
        OrdinaryDiffEq.ODEProblem(fn, u0, tspan)
    elseif jacobian === :fd
        # No Jacobian attached: the integrator differentiates the wrapped RHS
        # itself, which requires finite differences (Float64 workspace).
        if integrator isa KenCarp4
            integrator = KenCarp4(autodiff = AutoFiniteDiff())
        end
        OrdinaryDiffEq.ODEProblem(rhs!, u0, tspan)
    else
        error("jacobian must be :analytic, :approx, or :fd, got $(repr(jacobian))")
    end

    save_times = if saveat !== nothing
        saveat
    elseif !isempty(problem.output_times)
        problem.output_times
    else
        [tspan[2]]
    end

    step_kwargs = dt === nothing ? (;) : (dt = dt, adaptive = false)
    ode_sol = OrdinaryDiffEq.solve(
        ode_problem,
        integrator;
        abstol = abstol,
        reltol = reltol,
        saveat = save_times,
        maxiters = maxiters,
        step_kwargs...,
    )

    nt = length(ode_sol.t)
    T = Matrix{Float64}(undef, n, nt)
    ξ = Array{Float64,3}(undef, NC, n, nt)
    z = Matrix{Float64}(undef, n, nt)
    for m = 1:nt
        u = ode_sol.u[m]
        @inbounds for i = 1:n
            T[i, m] = u[i]
            z[i, m] = use_ale ? (u[zoff+i] + u[zoff+i+1]) / 2 : zc0[i]
            for j = 1:NC
                ξ[j, i, m] = u[n+(j-1)*n+i]
            end
        end
    end
    retcode =
        OrdinaryDiffEq.SciMLBase.successful_retcode(ode_sol) ? :Success :
        Symbol(ode_sol.retcode)
    return MMSSolution(retcode, collect(ode_sol.t), T, ξ, z)
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
