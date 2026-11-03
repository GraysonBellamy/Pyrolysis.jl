# ============================================================================
# Blowing: how gas leaving through a face changes its convective coefficient.
#
# The model types and their constants live in Materials/blowing_models.jl.
# This file evaluates them at a boundary face: the surface Newton builds a
# `BlowingFace` from the outward species fluxes at every iterate, and each
# `ConvectiveBC` asks its model for the ratio χ = h/h₀.
# ============================================================================

using ForwardDiff
using ..Materials: R_GAS, P_REF, evaluate_param
using ..Properties: chapman_enskog_diffusivity

"""
	BlowingFace{NC,T,M}

What a blowing model sees at a boundary face: the outward mass flux of every
component `ṁ` [kg/(m²·s), positive = leaving, zero for condensed components],
the surface-temperature iterate `T_s` [K] and the `material` whose gas
components supply the species properties. Built by the surface Newton at
every iterate; `T` is `Float64` in the residual and `ForwardDiff.Dual` under
AD.
"""
struct BlowingFace{NC,T<:Real,M<:Material,X}
    ṁ::NTuple{NC,T}
    T_s::T
    material::M
    surface_mole_fractions::X
end

@inline function BlowingFace(
    ṁ::NTuple{NC,<:Real},
    T_s::Real,
    material::Material{NC};
    surface_mole_fractions = nothing,
) where {NC}
    Tu = promote_type(eltype(ṁ), typeof(T_s))
    return BlowingFace{NC,Tu,typeof(material),typeof(surface_mole_fractions)}(
        ntuple(j -> Tu(ṁ[j]), Val(NC)),
        Tu(T_s),
        material,
        surface_mole_fractions,
    )
end

# ----------------------------------------------------------------------------
# The film factor
# ----------------------------------------------------------------------------

"""
	blowing_factor(φ) -> Real

Film (Couette, Spalding) blowing factor `φ/(e^φ − 1)`: the ratio of the heat
exchange across a still gas layer through which gas leaves at blowing
parameter `φ` to the exchange without the flow.

Smooth for all real `φ`: `f(0) = 1` (no flow), `f(φ) → φ·e^(−φ) → 0` as
`φ → ∞` (exchange extinguished by strong outgassing), `f(φ) > 1` for
`φ < 0` (suction enhancement). Three numerically stable branches:

- `|φ| < 1e-5`: 4th-order Taylor `1 − φ/2 + φ²/12` (the naive form is
  0/0 at the origin; values and first derivatives match the exact form
  to ~1e-16 at the branch seams).
- `φ > 0`: `φ·e^(−φ)/(−expm1(−φ))` — the naive `φ/expm1(φ)` overflows to
  `Inf` near `φ ≈ 710` and poisons ForwardDiff partials with `NaN`; this
  form underflows cleanly to zero with clean Duals.
- `φ < 0`: `φ/expm1(φ)` — `expm1(φ) → −1`, no overflow possible.

AD-clean: generic over `Real`, branches on the primal value only through
ordinary comparisons, so `ForwardDiff.Dual` inputs propagate exactly.
"""
@inline function blowing_factor(φ::Real)
    if abs(φ) < 1e-5
        return 1 - φ / 2 + φ^2 / 12
    elseif φ > 0
        return φ * exp(-φ) / (-expm1(-φ))
    else
        return φ / expm1(φ)
    end
end

# Net blowing conductance ṁ″·c_p,g = Σ_g ṁ_g · cp_g(T_s) [W/(m²·K)] at a
# boundary face, summed over gas components with the SIGNED outward mass flux
# (outflow > 0 pushes the boundary layer off — blowing; net inflow gives a
# negative conductance — suction — which `blowing_factor` maps to an
# enhancement). cp is evaluated at the current surface-temperature iterate;
# surface evaluation is a frozen-property thermal-film approximation. It
# keeps the residual a pure function of state. Compile-time unrolling removes
# the gas predicate from the generated evaluation.
@generated function _blowing_conductance(
    material::Material{NC},
    ṁ::NTuple{NC,<:Real},
    T_s::Real,
) where {NC}
    exprs = Expr[]
    push!(exprs, :(Tu = promote_type(eltype(ṁ), typeof(T_s))))
    push!(exprs, :(mcp = zero(Tu)))
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                mcp += ṁ[$j] * material.components[$j].heat_capacity(T_s)
            end
        end)
    end
    push!(exprs, :(return mcp))
    return Expr(:block, exprs...)
end

# Signed net outward gas flux through the face [kg/(m²·s)].
@generated function _gas_outflow(material::Material{NC}, ṁ::NTuple{NC,<:Real}) where {NC}
    exprs = Expr[]
    push!(exprs, :(tot = zero(eltype(ṁ))))
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                tot += ṁ[$j]
            end
        end)
    end
    push!(exprs, :(return tot))
    return Expr(:block, exprs...)
end

# Per-component gas properties at T, as one homogeneous tuple so the mixing
# loops below are type-stable: (gas?, cp, k, M, name). Condensed components get
# zeros and a unit molar mass they never divide by.
@generated function _gas_props(material::Material{NC}, T::Real) where {NC}
    exprs = Expr[]
    push!(exprs, :(Tz = zero(T)))
    items = Expr[]
    for j = 1:NC
        push!(
            items,
            quote
                @inbounds (
                    c = material.components[$j];
                    is_gas(c) ?
                    (
                        gas = true,
                        cp = c.heat_capacity(T) + Tz,
                        k = c.thermal_conductivity(T) + Tz,
                        M = molar_mass(c),
                        name = c.name,
                    ) : (gas = false, cp = Tz, k = Tz, M = 1.0, name = c.name)
                )
            end,
        )
    end
    push!(exprs, Expr(:return, Expr(:tuple, items...)))
    return Expr(:block, exprs...)
end

# ----------------------------------------------------------------------------
# Leaving-mixture properties
# ----------------------------------------------------------------------------

# Heat capacity and conductivity of the gas leaving the face. Lumped volatiles
# state them directly.
@inline function _leaving_mixture(vol::GasSpecies, face::BlowingFace)
    T = gas_property_temperature(vol, face.T_s)
    return vol.heat_capacity(T), vol.thermal_conductivity(T)
end

# Per-species: the components with outward flux, heat capacity weighted by
# mass fraction, conductivity by Wilke's rule at the leaving mole fractions
# with each species' viscosity from Eucken's relation (components carry none).
function _leaving_mixture(::Nothing, face::BlowingFace{NC}) where {NC}
    props = _gas_props(face.material, face.T_s)
    ṁ = face.ṁ
    Tu = eltype(ṁ)
    tot = zero(Tu)
    cp = zero(Tu)
    @inbounds for j = 1:NC
        if props[j].gas && ForwardDiff.value(ṁ[j]) > 0
            tot += ṁ[j]
            cp += ṁ[j] * props[j].cp
        end
    end
    cp /= tot
    k = zero(Tu)
    @inbounds for i = 1:NC
        (props[i].gas && ForwardDiff.value(ṁ[i]) > 0) || continue
        xi = ṁ[i] / props[i].M                 # ∝ mole fraction; the scale cancels
        μi = eucken_viscosity(props[i].k, props[i].cp, props[i].M)
        den = zero(Tu)
        for j = 1:NC
            (props[j].gas && ForwardDiff.value(ṁ[j]) > 0) || continue
            xj = ṁ[j] / props[j].M
            μj = eucken_viscosity(props[j].k, props[j].cp, props[j].M)
            Φ =
                (1 + sqrt(μi / μj) * (props[j].M / props[i].M)^0.25)^2 /
                sqrt(8 * (1 + props[i].M / props[j].M))
            den += xj * Φ
        end
        k += xi * props[i].k / den
    end
    return cp, k
end

# Blowing parameter of the property-ratio form, φ = ṁ″ c_p,∞ (Γ_∞/Γ_g)^n / h₀.
# Under net inflow the gas arriving is the free stream itself (ratio 1).
@inline function _property_ratio_phi(
    volatiles,
    face::BlowingFace,
    fs::GasSpecies,
    ṁ,
    cp∞,
    h,
    n_Γ,
)
    ForwardDiff.value(ṁ) > 0 || return ṁ * cp∞ / h
    cp_g, k_g = _leaving_mixture(volatiles, face)
    Γ∞ = fs.thermal_conductivity(gas_property_temperature(fs, face.T_s)) / cp∞
    Γg = k_g / cp_g
    return ṁ * cp∞ * (Γ∞ / Γg)^n_Γ / h
end

# ----------------------------------------------------------------------------
# Model evaluation
# ----------------------------------------------------------------------------

"""
	blowing_state(model, face, h, T_inf) -> (; χ, φ, ψ, Y, Λ, states)

Evaluate a blowing model at a `BlowingFace` against the bare coefficient `h`
and the ambient temperature `T_inf`: the ratio `χ = h_eff/h` the convective
boundary applies, the blowing parameter `φ`, and for `PoolingFilm` the pooling
factor `ψ`, the face mass fraction `Y` of the volatiles, the face buoyancy `Λ`
and the number of self-consistent `states` (1, or 2 where the mean of two is
taken). The other models report `ψ = 1`, `Y = 0`, `Λ = 1`, `states = 1`.
All entries are AD-clean.
"""
function blowing_state end

@inline function blowing_state(::NoBlowing, face::BlowingFace, h, T_inf)
    o = one(face.T_s)
    return (χ = o, φ = zero(o), ψ = o, Y = zero(o), Λ = o, states = 1)
end

@inline function blowing_state(m::CouetteFilm, face::BlowingFace, h, T_inf)
    mcp = _blowing_conductance(face.material, face.ṁ, face.T_s)
    φ = 2 * m.λ * mcp / h
    χ = blowing_factor(φ)
    o = one(χ)
    return (χ = χ, φ = φ, ψ = o, Y = zero(o), Λ = o, states = 1)
end

@inline function blowing_state(m::FunctionBlowing, face::BlowingFace, h, T_inf)
    χ = m.f(face, h, T_inf)
    o = one(χ)
    return (χ = χ, φ = zero(o), ψ = o, Y = zero(o), Λ = o, states = 1)
end

function blowing_state(m::PropertyRatioFilm, face::BlowingFace, h, T_inf)
    _empirical_signed_domain(m, face)
    fs = m.free_stream
    ṁ = _gas_outflow(face.material, face.ṁ)
    cp∞ = fs.heat_capacity(gas_property_temperature(fs, face.T_s))
    φ = _property_ratio_phi(m.volatiles, face, fs, ṁ, cp∞, h, m.n_Γ)
    χ = blowing_factor(φ)
    o = one(χ)
    return (χ = χ, φ = φ, ψ = o, Y = zero(o), Λ = o, states = 1)
end

"""
	blowing_ratio(model, face, h, T_inf) -> χ

The ratio `h_eff/h` of [`blowing_state`](@ref).
"""
@inline blowing_ratio(model::AbstractBlowingModel, face::BlowingFace, h, T_inf) =
    blowing_state(model, face, h, T_inf).χ

# ----------------------------------------------------------------------------
# Pooling model
# ----------------------------------------------------------------------------

# Five-point Gauss–Legendre mean of a property over [T_a, T_b]: exact for the
# polynomial heat capacities the gas data carry.
const _GL5_X =
    (-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640)
const _GL5_W = (
    0.2369268850561891,
    0.4786286704993665,
    0.5688888888888889,
    0.4786286704993665,
    0.2369268850561891,
)
@inline function _mean_property(p, T_a, T_b)
    Tm = (T_a + T_b) / 2
    dT = (T_b - T_a) / 2
    s = zero(Tm)
    @inbounds for i = 1:5
        s += _GL5_W[i] * p(Tm + dT * _GL5_X[i])
    end
    return s / 2
end

# The Fuller diffusion volume of component `name`, NaN without an entry: a
# search over the names, since a field read by a runtime name boxes its value.
@inline function _diffusion_volume(V::NamedTuple{names}, name::Symbol) where {names}
    vals = values(V)
    v = NaN
    for i in eachindex(names)
        names[i] === name && (v = vals[i])
    end
    return v
end

# Binary diffusivity of a gas component into the free stream [m²/s].
@inline function _component_diffusivity(m::PoolingFilm, name::Symbol, M, T)
    V_a = _diffusion_volume(m.diffusion_volume, name)
    if !isnan(V_a) && !isnan(m.free_stream.diffusion_volume)
        return fuller_diffusivity(
            T,
            M,
            V_a,
            m.free_stream.molar_mass,
            m.free_stream.diffusion_volume,
        )
    end
    return chapman_enskog_diffusivity(T, M)
end

@inline function _lumped_diffusivity(vol::GasSpecies, fs::GasSpecies, T)
    vol.diffusivity === nothing || return vol.diffusivity(gas_property_temperature(vol, T))
    if !isnan(vol.diffusion_volume) && !isnan(fs.diffusion_volume)
        return fuller_diffusivity(
            T,
            vol.molar_mass,
            vol.diffusion_volume,
            fs.molar_mass,
            fs.diffusion_volume,
        )
    end
    return chapman_enskog_diffusivity(T, vol.molar_mass)
end

# The species the face balance tracks: (ω, 1/M, a) with ω the share of the
# outflow, and a = ṁ″/(g s^{n_D}) the film parameter without the pooling
# factor. One entry for lumped volatiles; one per component otherwise, with
# ω = 0 for components that are condensed, free-stream gas or not leaving.
@inline function _face_species(m::PoolingFilm, vol::GasSpecies, face::BlowingFace, ṁ, g, α∞)
    s = _lumped_diffusivity(vol, m.free_stream, face.T_s) / α∞
    return ((ω = one(ṁ), Minv = 1 / vol.molar_mass, a = ṁ / (g * s^m.n_D)),)
end

function _face_species(
    m::PoolingFilm,
    ::Nothing,
    face::BlowingFace{NC},
    ṁ,
    g,
    α∞,
) where {NC}
    props = _gas_props(face.material, face.T_s)
    fs_name = m.free_stream.name
    T_s = face.T_s
    return ntuple(Val(NC)) do j
        p = props[j]
        blown = p.gas && p.name !== fs_name && ForwardDiff.value(face.ṁ[j]) > 0
        if blown
            s = _component_diffusivity(m, p.name, p.M, T_s) / α∞
            (ω = face.ṁ[j] / ṁ, Minv = 1 / p.M, a = ṁ / (g * s^m.n_D))
        else
            (ω = zero(ṁ), Minv = 1 / p.M, a = zero(ṁ))
        end
    end
end

# G(ψ) = Ψ(Λ(ψ)) − ψ with its derivative, the face fraction and the buoyancy
# at ψ. Face fractions Y_i = ω_i (1 − e^{−a_i/ψ}); M_face from them and the
# free stream; Λ = [1 − (M_face/T_s)(T_∞/M_∞)]/(1 − T_∞/T_s); Ψ the logistic
# switch between the floor and 1.
@inline function _pool_residual(ψ, species, Minv∞, T_s, T_inf, ψ0, Λc, Λw)
    S = Minv∞ + zero(ψ)
    dS = zero(S)
    Ysum = zero(S)
    for sp in species
        e = exp(-sp.a / ψ)
        Y = sp.ω * (1 - e)
        dY = -sp.ω * e * sp.a / ψ^2
        S += Y * (sp.Minv - Minv∞)
        dS += dY * (sp.Minv - Minv∞)
        Ysum += Y
    end
    M_face = 1 / S
    dM = -dS / S^2
    d = 1 - T_inf / T_s
    c = T_inf * Minv∞ / T_s
    Λ = (1 - c * M_face) / d
    dΛ = -c * dM / d
    z = clamp((Λ - Λc) / Λw, -60.0, 60.0)
    σ = 1 / (1 + exp(-z))
    Ψ = ψ0 + (1 - ψ0) * σ
    dΨ = (1 - ψ0) * σ * (1 - σ) / Λw * dΛ
    return Ψ - ψ, dΨ - 1, Ysum, Λ
end

# Newton with a bisection safeguard inside a sign-change bracket [a, b].
@inline function _refine_root(G, a, b, ga)
    lo, hi, glo = a, b, ga
    x = (a + b) / 2
    for _ = 1:80
        g, dg = G(x)
        abs(g) < 1e-14 && break
        if sign(g) == sign(glo)
            lo, glo = x, g
        else
            hi = x
        end
        xn = x - g / dg
        (lo < xn < hi) || (xn = (lo + hi) / 2)
        abs(xn - x) < 1e-15 && (x = xn; break)
        x = xn
    end
    return x
end

# Primal roots of G on [ψ0, 1]: the lowest and highest self-consistent
# pooling factors (equal when there is one), and their count. G(ψ0) > 0 and
# G(1) < 0 by construction, so at least one root exists.
function _pool_roots(species, Minv∞, T_s, T_inf, ψ0, Λc, Λw)
    v = ForwardDiff.value
    sp = map(s -> (ω = v(s.ω), Minv = s.Minv, a = v(s.a)), species)
    Tsv, Tiv = v(T_s), v(T_inf)
    G = ψ -> begin
        g, dg, _, _ = _pool_residual(ψ, sp, Minv∞, Tsv, Tiv, ψ0, Λc, Λw)
        (g, dg)
    end
    N = 16
    ψprev = ψ0
    gprev, _ = G(ψprev)
    lo = hi = NaN
    n = 0
    for k = 1:N
        ψk = ψ0 + (1 - ψ0) * k / N
        gk, _ = G(ψk)
        if gk == 0
            r = ψk
        elseif gprev * gk < 0
            r = _refine_root(G, ψprev, ψk, gprev)
        else
            ψprev, gprev = ψk, gk
            continue
        end
        n == 0 && (lo = r)
        hi = r
        n += 1
        ψprev, gprev = ψk, gk
    end
    n == 0 && (lo =  hi = 1.0 ; n = 1)
    return lo, hi, n
end

function blowing_state(m::PoolingFilm, face::BlowingFace, h, T_inf)
    _empirical_signed_domain(m, face)
    names=map(c->c.name, face.material.components)
    fs = m.free_stream
    T_s = face.T_s
    T_p = gas_property_temperature(fs, T_s)      # property evaluations stay inside the fit range
    ṁ = _gas_outflow(face.material, face.ṁ)
    cp∞ = fs.heat_capacity(T_p)
    φ = _property_ratio_phi(m.volatiles, face, fs, ṁ, cp∞, h, m.n_Γ)
    o = one(φ)
    d = 1 - T_inf / T_s
    # No outflow, or a face no hotter than the free stream: nothing pools.
    if ForwardDiff.value(ṁ) <= 0 || ForwardDiff.value(d) <= 1e-6
        return (χ = blowing_factor(φ), φ = φ, ψ = o, Y = zero(o), Λ = o, states = 1)
    end
    if face.surface_mole_fractions !== nothing
        x=face.surface_mole_fractions
        masses=map(molar_mass, face.material.components)
        Mface=sum(j->x[j]*masses[j], eachindex(x))
        Λ=(1-Mface*T_inf/(fs.molar_mass*T_s))/d
        z=clamp((Λ-m.Λc)/m.Λw, -60.0, 60.0)
        ψ=m.ψ0+(1-m.ψ0)/(1+exp(-z))
        Y=sum(j->names[j]===fs.name ? zero(φ) : x[j]*masses[j]/Mface, eachindex(x))
        return (χ = ψ*blowing_factor(φ/ψ), φ = φ, ψ = ψ, Y = Y, Λ = Λ, states = 1)
    end
    # Without a solved mass film, this empirical closure assumes zero free-
    # stream species flux. Product shares then sum to one by species balance.
    for j in eachindex(face.ṁ)
        if names[j]===fs.name && ForwardDiff.value(face.ṁ[j])!=0
            throw(
                ArgumentError(
                    "PoolingFilm without a solved surface composition requires zero free-stream species flux",
                ),
            )
        end
    end
    # Gas exchange without outflow, and the film parameter of each species.
    g = h / _mean_property(fs.heat_capacity, gas_property_temperature(fs, T_inf), T_p)
    ρ∞ = P_REF * fs.molar_mass / (R_GAS * T_p)
    α∞ = fs.thermal_conductivity(T_p) / (ρ∞ * cp∞)
    species = _face_species(m, m.volatiles, face, ṁ, g, α∞)
    Minv∞ = 1 / fs.molar_mass
    # Self-consistent ψ on the primal values, then one implicit-function step in
    # the working scalar type so AD partials ride through the root.
    lo, hi, n = _pool_roots(species, Minv∞, T_s, T_inf, m.ψ0, m.Λc, m.Λw)
    lift(r) = begin
        g_r, dg_r, _, _ =
            _pool_residual(r + zero(φ), species, Minv∞, T_s, T_inf, m.ψ0, m.Λc, m.Λw)
        r - g_r / dg_r
    end
    ψ_hi = lift(hi)
    χ_hi = ψ_hi * blowing_factor(φ / ψ_hi)
    if n > 1
        ψ_lo = lift(lo)
        χ_lo = ψ_lo * blowing_factor(φ / ψ_lo)
        χ = (χ_lo + χ_hi) / 2
        ψ = (ψ_lo + ψ_hi) / 2
    else
        χ = χ_hi
        ψ = ψ_hi
    end
    _, _, Y, Λ = _pool_residual(ψ_hi, species, Minv∞, T_s, T_inf, m.ψ0, m.Λc, m.Λw)
    return (χ = χ, φ = φ, ψ = ψ, Y = Y, Λ = Λ, states = n)
end

# ----------------------------------------------------------------------------
# Diagnostics: the blowing state of the first convective member of a thermal BC
# ----------------------------------------------------------------------------

_face_blowing_state(::AbstractThermalBC, face::BlowingFace, t) = nothing

function _face_blowing_state(bc::ConvectiveBC, face::BlowingFace, t)
    h = evaluate_param(bc.h, t)
    T_inf = evaluate_param(bc.T_inf, t)
    ForwardDiff.value(h) > 0 || return nothing
    return blowing_state(bc.blowing, face, h, T_inf)
end

@generated function _face_blowing_state(
    bc::CombinedThermalBC{BCs},
    face::BlowingFace,
    t,
) where {BCs}
    exprs = Expr[]
    for i = 1:fieldcount(BCs)
        push!(exprs, quote
            s = _face_blowing_state(bc.bcs[$i], face, t)
            s === nothing || return s
        end)
    end
    push!(exprs, :(return nothing))
    return Expr(:block, exprs...)
end

# Write the converged blowing diagnostics onto a BoundaryState (Float64).
function _record_blowing!(state::BoundaryState, bc::AbstractThermalBC, face::BlowingFace, t)
    v = ForwardDiff.value
    state.blowing_mcp = v(_blowing_conductance(face.material, face.ṁ, face.T_s))
    s = _face_blowing_state(bc, face, t)
    if s === nothing
        state.blowing_ratio = 1.0
        state.pooling_factor = 1.0
        state.face_fraction = 0.0
        state.face_buoyancy = 1.0
    else
        state.blowing_ratio = v(s.χ)
        state.pooling_factor = v(s.ψ)
        state.face_fraction = v(s.Y)
        state.face_buoyancy = v(s.Λ)
    end
    return nothing
end

# Empirical property-ratio/pooling data support product outflow and pure
# free-stream suction. Mixed counterflow has no derived empirical extension;
# the signed Couette film supplies a separately derived closure for that regime.
function _empirical_signed_domain(model, face)
    names=map(c->c.name, face.material.components)
    incoming=false
    outgoing=false
    for j in eachindex(face.ṁ)
        v=ForwardDiff.value(face.ṁ[j])
        incoming |= v<0
        outgoing |= v>0
        if v<0 && names[j]!==model.free_stream.name
            throw(
                ArgumentError(
                    "$(typeof(model).name.name) supports suction of its free_stream species only",
                ),
            )
        end
    end
    incoming &&
        outgoing &&
        throw(
            ArgumentError(
                "$(typeof(model).name.name) has no qualified mixed-ingress/outflow correlation",
            ),
        )
    return nothing
end
