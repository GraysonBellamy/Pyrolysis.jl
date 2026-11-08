# =============================================================================
# Exact solution: the gray two-stream rows of a slab
# =============================================================================
#
# The solver's radiative transfer rows, solved exactly for given cell optical
# thicknesses and emissive powers. Cells 1..N run from the exposed face (face
# 1) to the back (face N + 1); the fluxes are hemispherical, per unit face
# area [W/m²]: q⁻ toward the back, q⁺ toward the exposed face. The
# hemispherical form of the transfer equation (Modest & Mazumder, Radiative
# Heat Transfer, 4th ed., §14.3) with the diffusivity D in place of the
# two-flux factor 2, integrated exactly across a cell of uniform emissive power
# Bᵢ and optical thickness τᵢ, gives with tᵢ = e^{−Dτᵢ}
#
#   q⁻ᵢ₊₁ = Bᵢ + (q⁻ᵢ − Bᵢ) tᵢ,        q⁺ᵢ = Bᵢ + (q⁺ᵢ₊₁ − Bᵢ) tᵢ.
#
# For a linear B the continuum limit of these rows carries −(2/(Dκ)) dB/dx, the
# Rosseland flux when D = 3/2; per cell drop g of B the discrete rows carry
# g·coth(Dτ/2) (substitute q± = B + const into the rows). Beams (a collimated
# heater, or the directions of a cone) travel along their direction cosine μ:
# I_{h,i+1} = I_{h,i} e^{−τᵢ/μ}.
#
# Faces are gray interfaces. The exposed face has the transmissivity ε both
# ways, so with Q the diffuse flux arriving from outside
#
#   q⁻₁ = ε Q + (1 − ε) q⁺₁,
#
# and each beam enters with its given flux. The back has the transmissivity
# ε_b (zero is a mirror) and the diffuse flux Q_b arriving from outside:
#
#   q⁺_b = ε_b Q_b + (1 − ε_b)(q⁻_b + Σ_h I_{h,b}),
#
# so the back reflects 1 − ε_b of every arriving beam into q⁺, where it travels
# on as diffuse radiation, and ε_b of it leaves. At a mirror the beams may
# instead reflect specularly into an upward family that returns to the exposed
# face, where 1 − ε of it reflects down again, closed as a geometric series.
#
# The thermal rows are linear in the fluxes. With q⁻₁ = a as a parameter one
# downward and one upward pass give q⁺₁ = A a + b, the exposed row gives a, and
# two more passes write the values; the determinant 1 − (1 − ε)(1 − ε_b)Πtᵢ²
# is positive for ε > 0. Cells before `first` are not part of the medium:
# their faces hold zero fluxes and face `first` is the exposed face.

"""
A beam family: direction cosine `μ` in the medium and the flux `q0` [W/m², per
unit face area] entering through the exposed face.
"""
struct TwoStreamBeam
    μ::Float64
    q0::Float64
end

"""
    two_stream_cone_beams(q, θ1, θ2; directions = 3, a = 1.0)

The beam families of a heater of flux `q` whose intensity is uniform between
the polar angles `θ1 < θ2` [rad], entering at the absorptivity `a`:
Gauss–Legendre nodes in μ over [cos θ2, cos θ1], each family carrying the share
of the flux that its weight times μ gives. The shares sum to one.
"""
function two_stream_cone_beams(q, θ1, θ2; directions::Int = 3, a = 1.0)
    μ1, μ2 = cos(θ1), cos(θ2)
    x, w = _gauss_legendre(directions)
    mid, half = (μ1 + μ2) / 2, (μ1 - μ2) / 2
    total = (μ1^2 - μ2^2) / 2
    return [
        TwoStreamBeam(mid + half * x[k], a * q * half * w[k] * (mid + half * x[k]) / total)
        for k in eachindex(x)
    ]
end

"""
Storage of a two-stream solution of `N` cells and `nb` beam families:
`t[i] = e^{−Dτᵢ}`, the face fluxes `qm` (q⁻) and `qp` (q⁺), and the beams `Id`
(downward) and `Iu` (upward, a specular mirror only), face by family.
"""
struct TwoStreamWork
    t::Vector{Float64}
    qm::Vector{Float64}
    qp::Vector{Float64}
    Id::Matrix{Float64}
    Iu::Matrix{Float64}
end
TwoStreamWork(N::Int, nb::Int = 0) =
    TwoStreamWork(zeros(N), zeros(N + 1), zeros(N + 1), zeros(N + 1, nb), zeros(N + 1, nb))

"""
    two_stream_sweep!(w, τ, B; D = 1.5, ε = 1.0, Q = 0.0, ε_b = 0.0, Q_b = 0.0,
                      beams = (), beam_back = :diffuse, first = 1) -> w

The fluxes that satisfy every row, written into `w` without allocating.
`beam_back = :diffuse` reflects the beams at the back into q⁺; `:specular`
(a mirror, `ε_b = 0`) into an upward beam family.
"""
function two_stream_sweep!(
    w::TwoStreamWork,
    τ,
    B;
    D = 1.5,
    ε = 1.0,
    Q = 0.0,
    ε_b = 0.0,
    Q_b = 0.0,
    beams = (),
    beam_back::Symbol = :diffuse,
    first::Int = 1,
)
    N = length(τ)
    s = first
    t, qm, qp, Id, Iu = w.t, w.qm, w.qp, w.Id, w.Iu
    specular = beam_back === :specular
    specular ||
        beam_back === :diffuse ||
        throw(ArgumentError("beam_back must be :diffuse or :specular"))
    specular &&
        ε_b != 0 &&
        throw(ArgumentError("specular beam reflection needs a mirror back (ε_b = 0)"))
    @inbounds for j = 1:(s-1)
        qm[j] = 0.0
        qp[j] = 0.0
        for h in axes(Id, 2)
            Id[j, h] = 0.0
            Iu[j, h] = 0.0
        end
    end
    @inbounds for i = s:N
        t[i] = exp(-D * τ[i])
    end

    # Beams, down from the exposed face; what reaches the back.
    arriving = 0.0
    @inbounds for (h, beam) in enumerate(beams)
        x = beam.q0
        if specular
            Tr = 1.0
            for i = s:N
                Tr *= exp(-τ[i] / beam.μ)
            end
            x /= 1 - (1 - ε) * Tr^2
        end
        Id[s, h] = x
        for i = s:N
            Id[i+1, h] = Id[i, h] * exp(-τ[i] / beam.μ)
        end
        if specular
            Iu[N+1, h] = Id[N+1, h]
            for i = N:-1:s
                Iu[i, h] = Iu[i+1, h] * exp(-τ[i] / beam.μ)
            end
        else
            for j = s:(N+1)
                Iu[j, h] = 0.0
            end
            arriving += Id[N+1, h]
        end
    end

    # Thermal field: every face affine in a = q⁻ at the exposed face.
    A, b = 1.0, 0.0
    @inbounds for i = s:N
        A, b = A * t[i], B[i] + (b - B[i]) * t[i]
    end
    ρ_b = 1 - ε_b
    A, b = ρ_b * A, ε_b * Q_b + ρ_b * (b + arriving)
    @inbounds for i = N:-1:s
        A, b = A * t[i], B[i] + (b - B[i]) * t[i]
    end
    qm[s] = (ε * Q + (1 - ε) * b) / (1 - (1 - ε) * A)
    @inbounds for i = s:N
        qm[i+1] = B[i] + (qm[i] - B[i]) * t[i]
    end
    qp[N+1] = ε_b * Q_b + ρ_b * (qm[N+1] + arriving)
    @inbounds for i = N:-1:s
        qp[i] = B[i] + (qp[i+1] - B[i]) * t[i]
    end
    return w
end

"""
    two_stream_sources(qm, qp, Id, Iu; first = 1) -> S

Cell sources [W/m², per unit face area] of the face fluxes, zero before the
cell `first`:
`Sᵢ = (q⁻ᵢ − q⁻ᵢ₊₁) + (q⁺ᵢ₊₁ − q⁺ᵢ) + Σ_h (I↓_{h,i} − I↓_{h,i+1} + I↑_{h,i+1} − I↑_{h,i})`.
"""
function two_stream_sources(qm, qp, Id, Iu; first::Int = 1)
    N = length(qm) - 1
    return [
        i < first ? 0.0 :
        (qm[i] - qm[i+1]) +
        (qp[i+1] - qp[i]) +
        sum(Id[i, h] - Id[i+1, h] + Iu[i+1, h] - Iu[i, h] for h in axes(Id, 2); init = 0.0)
        for i = 1:N
    ]
end

"""
    two_stream_solve(τ, B; D, ε, Q, ε_b, Q_b, beams, beam_back, first)
        -> (; qm, qp, Id, Iu, S, q_out, q_enter, q_back)

The fluxes of `two_stream_sweep!` with the cell sources `S`, the flux leaving
through the exposed face (`ε (q⁺ + Σ I↑)` there), the flux entering through
it (`ε Q + Σ_h q0_h`) and the net flux leaving through the back
(`q⁻ + Σ I↓ − q⁺ − Σ I↑` there).
"""
function two_stream_solve(τ, B; beams = (), ε = 1.0, Q = 0.0, first::Int = 1, kwargs...)
    N, nb = length(τ), length(beams)
    w = two_stream_sweep!(TwoStreamWork(N, nb), τ, B; beams, ε, Q, first, kwargs...)
    qm, qp, Id, Iu = w.qm, w.qp, w.Id, w.Iu
    up(j) = sum(Iu[j, h] for h = 1:nb; init = 0.0)
    down(j) = sum(Id[j, h] for h = 1:nb; init = 0.0)
    return (;
        qm,
        qp,
        Id,
        Iu,
        S = two_stream_sources(qm, qp, Id, Iu; first),
        q_out = ε * (qp[first] + up(first)),
        q_enter = ε * Q + sum(b.q0 for b in beams; init = 0.0),
        q_back = qm[N+1] + down(N + 1) - qp[N+1] - up(N + 1),
    )
end

"""
    two_stream_residual(qm, qp, Id, Iu, τ, B; D, ε, Q, ε_b, Q_b, beams, beam_back, first)
        -> (; thermal, down, up)

Every row at trial fluxes, each as the flux minus what its row gives:
`thermal[1, j]` the q⁻ row of face j, `thermal[2, j]` its q⁺ row, `down[j, h]`
and `up[j, h]` the rows of beam h at face j. A face before `first`, and an
upward beam unless the mirror reflects specularly, holds its own value.
"""
function two_stream_residual(
    qm,
    qp,
    Id,
    Iu,
    τ,
    B;
    D = 1.5,
    ε = 1.0,
    Q = 0.0,
    ε_b = 0.0,
    Q_b = 0.0,
    beams = (),
    beam_back::Symbol = :diffuse,
    first::Int = 1,
)
    N, nb, s = length(τ), length(beams), first
    T = promote_type(eltype(qm), eltype(qp), eltype(Id), eltype(Iu), eltype(τ), eltype(B))
    thermal, down, up = zeros(T, 2, N + 1), zeros(T, N + 1, nb), zeros(T, N + 1, nb)
    specular = beam_back === :specular
    for j = 1:(s-1)
        thermal[1, j], thermal[2, j] = qm[j], qp[j]
        for h = 1:nb
            down[j, h], up[j, h] = Id[j, h], Iu[j, h]
        end
    end
    thermal[1, s] = qm[s] - (ε * Q + (1 - ε) * qp[s])
    for (h, beam) in enumerate(beams)
        down[s, h] = Id[s, h] - (beam.q0 + (specular ? (1 - ε) * Iu[s, h] : zero(T)))
    end
    for i = s:N
        ti = exp(-D * τ[i])
        thermal[1, i+1] = qm[i+1] - (B[i] + (qm[i] - B[i]) * ti)
        thermal[2, i] = qp[i] - (B[i] + (qp[i+1] - B[i]) * ti)
        for (h, beam) in enumerate(beams)
            th = exp(-τ[i] / beam.μ)
            down[i+1, h] = Id[i+1, h] - Id[i, h] * th
            up[i, h] = specular ? Iu[i, h] - Iu[i+1, h] * th : Iu[i, h]
        end
    end
    arriving = specular ? zero(T) : sum(Id[N+1, h] for h = 1:nb; init = zero(T))
    thermal[2, N+1] = qp[N+1] - (ε_b * Q_b + (1 - ε_b) * (qm[N+1] + arriving))
    for h = 1:nb
        up[N+1, h] = specular ? Iu[N+1, h] - Id[N+1, h] : Iu[N+1, h]
    end
    return (; thermal, down, up)
end

# -----------------------------------------------------------------------------
# Closed forms
# -----------------------------------------------------------------------------

"""
    two_stream_isothermal_slab(; B, x, ε, Q = 0.0, ε_b, Q_b = 0.0)

The rows over a slab of uniform emissive power `B` and total `Dτ_L = x`, with
the diffuse inputs `Q` (exposed face, transmissivity `ε`) and `Q_b` (back,
transmissivity `ε_b`, zero for a mirror). With y = q⁻₁ − B and z = q⁺_b − B,
t = e^{−x}, the exposed and back rows read y = ε(Q − B) + (1 − ε) t z and
z = ε_b(Q_b − B) + (1 − ε_b) t y. Returns the four face fluxes `qm_e`, `qp_e`,
`qm_b`, `qp_b`, the net intake `q⁻₁ − q⁺₁`, the flux leaving through the exposed
face `ε q⁺₁` and the net flux leaving through the back `q⁻_b − q⁺_b`.
"""
function two_stream_isothermal_slab(; B, x, ε, Q = 0.0, ε_b, Q_b = 0.0)
    t = exp(-x)
    y = (ε * (Q - B) + (1 - ε) * ε_b * t * (Q_b - B)) / (1 - (1 - ε) * (1 - ε_b) * t^2)
    z = ε_b * (Q_b - B) + (1 - ε_b) * t * y
    qm_e, qp_e, qm_b, qp_b = B + y, B + z * t, B + y * t, B + z
    return (;
        qm_e,
        qp_e,
        qm_b,
        qp_b,
        intake = y - z * t,
        q_out = ε * qp_e,
        q_back = y * t - z,
    )
end

"""
    two_stream_slab_emission(B, x, ε, ε_b)

The flux an isothermal slab of emissive power `B` and total `Dτ_L = x` emits
through its exposed face (transmissivity `ε`) into surroundings at 0 K, behind
it a back of transmissivity `ε_b` facing 0 K (zero: a mirror):
`εB(1 − t)(1 + (1 − ε_b)t) / (1 − (1 − ε)(1 − ε_b)t²)`, t = e^{−x}, with
`1 − t = −expm1(−x)` so a thin slab keeps its digits.
"""
function two_stream_slab_emission(B, x, ε, ε_b)
    t = exp(-x)
    return ε * B * (-expm1(-x)) * (1 + (1 - ε_b) * t) / (1 - (1 - ε) * (1 - ε_b) * t^2)
end

"""
    two_stream_linear_flux(g, x)

The net flux toward the back at every interior face when the emissive power
drops by `g` per cell of `Dτ = x`: `g·coth(x/2)`, evaluated as
`g(2 + expm1(−x))/(−expm1(−x))`. With q± = B + const the cell rows give
q⁻ − B = g/(1 − t) and q⁺ − B = −g t/(1 − t) at the face upstream of a cell.
"""
two_stream_linear_flux(g, x) = g * (2 + expm1(-x)) / (-expm1(-x))

"""
    two_stream_beam_slab(q0, τ_L, ε, ε_b; D = 1.5, μ = 1.0)

A beam of flux `q0` entering, along μ, a slab with no emission and no diffuse
input, of total optical thickness `τ_L`; the back (transmissivity `ε_b`, zero
for a mirror) reflects 1 − ε_b of it into q⁺. With t = e^{−Dτ_L} and
T_c = e^{−τ_L/μ}: q⁺₁ = t(1 − ε_b) q0 T_c / (1 − (1 − ε)(1 − ε_b)t²). Returns
the flux leaving through the exposed face `ε q⁺₁`, the flux leaving through
the back `ε_b (q⁻_b + q0 T_c)` and the absorbed remainder.
"""
function two_stream_beam_slab(q0, τ_L, ε, ε_b; D = 1.5, μ = 1.0)
    t, T_c = exp(-D * τ_L), exp(-τ_L / μ)
    qp_e = t * (1 - ε_b) * q0 * T_c / (1 - (1 - ε) * (1 - ε_b) * t^2)
    q_out = ε * qp_e
    q_back = ε_b * (t * (1 - ε) * qp_e + q0 * T_c)
    return (; q_out, q_back, absorbed = q0 - q_out - q_back)
end
