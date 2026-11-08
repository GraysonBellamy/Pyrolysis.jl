# =============================================================================
# Exact solution: gray radiative transfer in a slab of piecewise-constant cells
# =============================================================================
#
# The radiative transfer equation of a non-scattering gray medium of
# refractive index one (Modest & Mazumder, Radiative Heat Transfer, 4th ed.,
# §13.3, Eqs. 13.29–13.35), solved exactly for a slab of N cells numbered from
# the exposed face (face 1) to the back (face N + 1). Cell i has the optical
# thickness τᵢ = κᵢΔᵢ and a uniform emissive power Bᵢ = σTᵢ⁴ [W/m²]: the
# emission is the finite-volume staircase. Each direction μ = cos θ ∈ (0, 1] is
# independent, and across cell i the intensity is exactly
#
#   I_out = I_in e^{−τᵢ/μ} + (Bᵢ/π)(1 − e^{−τᵢ/μ}).
#
# Each face is a gray interface: transmissivity ε at the exposed face and ε_b
# at the back, the same both ways, reflectivity 1 − ε, the reflection specular
# or diffuse. Specular faces couple each direction to itself through a 2×2
# solve; a diffuse reflector couples the directions through one radiosity,
# closed by linear superposition. Outside each face the surroundings are
# diffuse at the emissive power B_amb (exposed) or B_b (back). The heater's
# flux q arrives either with uniform intensity between the polar angles θ₁ and
# θ₂ (diffuse: 0 and π/2), or as a collimated normal beam, which reflects
# specularly at a specular face and feeds the radiosity of a diffuse one.
#
# The angular integrals are composite Gauss–Legendre over sub-intervals of μ
# that break at cos θ₁ and cos θ₂, where the heater's intensity jumps; every
# integrand is then smooth on each sub-interval and the quadrature error is
# about 1e-14 at the default order. Fluxes are per unit face area: q is the
# net flux toward the back at each face, and cell i's source is qᵢ − qᵢ₊₁.
#
# Two closed forms need no angular quadrature. With black faces (ε = ε_b = 1),
# Modest & Mazumder's Eq. 13.30 integrated by parts (dE₃/dx = −E₂) gives the
# net flux at the optical depth τ as a sum over the jumps J_k of the emissive
# power at the depths τ_k, the surroundings standing as black half-spaces at
# B_amb and B_b:
#
#   q(τ) = 2q E₃(τ) − 2 Σ_k J_k E₃(|τ − τ_k|),
#
# a diffuse heater of flux q included. A heater of uniform intensity between θ₁
# and θ₂ crossing a cold, black-faced slab carries, with μ = cos θ and
# ∫₀^c μ e^{−τ/μ} dμ = c² E₃(τ/c),
#
#   q(τ) = 2q [μ₁² E₃(τ/μ₁) − μ₂² E₃(τ/μ₂)] / (μ₁² − μ₂²).

using SpecialFunctions: expint

"""
    gray_slab_black_faces(τ, B; q = 0.0, B_amb = 0.0, B_b = 0.0) -> q_faces

The net flux toward the back at each face of the slab with black faces, from
the jumps of the emissive power (closed form above); `q` is a diffuse heater.
"""
function gray_slab_black_faces(τ, B; q = 0.0, B_amb = 0.0, B_b = 0.0)
    N = length(τ)
    depth = [0.0; cumsum(τ)]
    jumps = [B[1] - B_amb; [B[i+1] - B[i] for i = 1:(N-1)]; B_b - B[N]]
    E3(x) = expint(3, x)
    return [
        2q * E3(depth[j]) - 2sum(jumps[k] * E3(abs(depth[j] - depth[k])) for k = 1:(N+1))
        for j = 1:(N+1)
    ]
end

"""
    gray_slab_cone_transmission(q, τ, θ1, θ2)

The flux of a heater `q` of uniform intensity between the polar angles
`θ1 < θ2` [rad] that crosses the optical depth `τ` of a cold slab with black
faces (closed form above); `θ1 = 0, θ2 = π/2` is the diffuse 2qE₃(τ).
"""
function gray_slab_cone_transmission(q, τ, θ1, θ2)
    μ1, μ2 = cos(θ1), cos(θ2)
    return 2q * (μ1^2 * expint(3, τ / μ1) - μ2^2 * expint(3, τ / μ2)) / (μ1^2 - μ2^2)
end

"""
    gray_slab_quadrature(; m = 24, cone = nothing) -> (μ, w)

Composite Gauss–Legendre nodes `μ ∈ (0, 1)` and weights, `m` per
sub-interval, with breakpoints at `cos θ₁` and `cos θ₂` when `cone = (θ₁, θ₂)`.
"""
function gray_slab_quadrature(; m::Int = 24, cone = nothing)
    breaks = [0.0, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5, 1.0]
    if cone !== nothing
        for c in (cos(cone[1]), cos(cone[2]))
            1e-12 < c < 1 - 1e-12 && push!(breaks, c)
        end
    end
    sort!(unique!(breaks))
    x, wx = _gauss_legendre(m)
    μ, w = Float64[], Float64[]
    for k = 1:(length(breaks)-1)
        a, b = breaks[k], breaks[k+1]
        append!(μ, (a + b) / 2 .+ (b - a) / 2 .* x)
        append!(w, (b - a) / 2 .* wx)
    end
    return μ, w
end

"""
    gray_slab_exact(τ, B; ε = 1.0, top_reflection = :specular, q = 0.0,
                    heater = :diffuse, B_amb = 0.0, ε_b = 1.0, B_b = 0.0,
                    back_reflection = :specular, m = 24)

Exact radiative transfer through the slab whose cells have the optical
thicknesses `τ` and the emissive powers `B` [W/m²], from the exposed face
down. `ε`, `ε_b`: the transmissivities of the exposed face and the back;
`top_reflection`, `back_reflection`: `:specular` or `:diffuse`. `heater`:
`:diffuse`, `:collimated` (a normal beam), or `(θ₁, θ₂)`, uniform intensity
between those polar angles [rad]. `B_amb`, `B_b`: the diffuse surroundings
outside the exposed face and the back.

Returns `(q, S, q_out, q_enter, q_back)`: the net flux toward the back at each
face, the cell sources `Sᵢ = qᵢ − qᵢ₊₁`, the flux leaving through the exposed
face, the flux entering through it (`ε (q + B_amb)`) and the net flux leaving
through the back (`q[end]`).
"""
function gray_slab_exact(
    τ,
    B;
    ε = 1.0,
    top_reflection::Symbol = :specular,
    q = 0.0,
    heater = :diffuse,
    B_amb = 0.0,
    ε_b = 1.0,
    B_b = 0.0,
    back_reflection::Symbol = :specular,
    m::Int = 24,
)
    N = length(τ)
    length(B) == N || throw(ArgumentError("τ and B need one value per cell"))
    for r in (top_reflection, back_reflection)
        r in (:specular, :diffuse) ||
            throw(ArgumentError("reflection must be :specular or :diffuse, got $r"))
    end
    cone =
        heater === :diffuse ? (0.0, π / 2) :
        heater === :collimated ? nothing : (Float64(heater[1]), Float64(heater[2]))
    μ, w = gray_slab_quadrature(; m, cone)
    K = length(μ)
    hemi(I) = 2π * sum(w[k] * μ[k] * I[k] for k = 1:K)

    # Specular and diffuse reflectivities of each face.
    Rt, ρt = top_reflection === :specular ? (1 - ε, 0.0) : (0.0, 1 - ε)
    Rb, ρb = back_reflection === :specular ? (1 - ε_b, 0.0) : (0.0, 1 - ε_b)

    # Intensity entering through each face from outside, per direction.
    μ1, μ2 = cone === nothing ? (1.0, 1.0) : (cos(cone[1]), cos(cone[2]))
    I_h = cone === nothing ? 0.0 : q / (π * (μ1^2 - μ2^2))
    St = [ε * (B_amb / π + (μ2 < μ[k] < μ1 ? I_h : 0.0)) for k = 1:K]
    Sb = ε_b * B_b / π

    # The collimated beam: down from the exposed face, specular reflections at
    # both faces summed as a geometric series.
    bd, bu = zeros(N + 1), zeros(N + 1)
    if cone === nothing
        Tr1 = exp(-sum(τ))
        bd[1] = ε * q / (1 - Rt * Rb * Tr1^2)
        for i = 1:N
            bd[i+1] = bd[i] * exp(-τ[i])
        end
        bu[N+1] = Rb * bd[N+1]
        for i = N:-1:1
            bu[i] = bu[i+1] * exp(-τ[i])
        end
    end

    # Every direction for the diffuse radiosities Jt, Jb [W/(m² sr)] added to
    # what enters the medium through the exposed face and the back.
    function sweep(Jt, Jb)
        Id, Iu = zeros(N + 1, K), zeros(N + 1, K)
        for k = 1:K
            Tr, Sdown, Sup = 1.0, 0.0, 0.0
            for i = 1:N
                t = exp(-τ[i] / μ[k])
                Tr *= t
                Sdown = Sdown * t - B[i] / π * expm1(-τ[i] / μ[k])
            end
            for i = N:-1:1
                Sup = Sup * exp(-τ[i] / μ[k]) - B[i] / π * expm1(-τ[i] / μ[k])
            end
            st, sb = St[k] + Jt, Sb + Jb
            det = 1 - Rt * Rb * Tr^2
            Id[1, k] = (st + Rt * Sup + Rt * Tr * (sb + Rb * Sdown)) / det
            Iu[N+1, k] = (sb + Rb * Sdown + Rb * Tr * (st + Rt * Sup)) / det
            for i = 1:N
                Id[i+1, k] = Id[i, k] * exp(-τ[i] / μ[k]) - B[i] / π * expm1(-τ[i] / μ[k])
            end
            for i = N:-1:1
                Iu[i, k] = Iu[i+1, k] * exp(-τ[i] / μ[k]) - B[i] / π * expm1(-τ[i] / μ[k])
            end
        end
        return Id, Iu
    end

    # The radiosities a sweep implies, the beam's arrivals included.
    implied(Id, Iu) = (
        ρt * (hemi(view(Iu, 1, :)) + bu[1]) / π,
        ρb * (hemi(view(Id, N + 1, :)) + bd[N+1]) / π,
    )
    Id, Iu = sweep(0.0, 0.0)
    if ρt > 0 || ρb > 0
        g0 = implied(Id, Iu)
        g1 = implied(sweep(1.0, 0.0)...)
        g2 = implied(sweep(0.0, 1.0)...)
        # J = g0 + M J, columns of M from the unit sweeps.
        a11, a21 = g1[1] - g0[1], g1[2] - g0[2]
        a12, a22 = g2[1] - g0[1], g2[2] - g0[2]
        d = (1 - a11) * (1 - a22) - a12 * a21
        Jt = ((1 - a22) * g0[1] + a12 * g0[2]) / d
        Jb = (a21 * g0[1] + (1 - a11) * g0[2]) / d
        Id, Iu = sweep(Jt, Jb)
    end

    qf = [hemi(view(Id, j, :)) - hemi(view(Iu, j, :)) + bd[j] - bu[j] for j = 1:(N+1)]
    S = [qf[i] - qf[i+1] for i = 1:N]
    q_out = ε * (hemi(view(Iu, 1, :)) + bu[1])
    return (; q = qf, S, q_out, q_enter = ε * (q + B_amb), q_back = qf[end])
end
