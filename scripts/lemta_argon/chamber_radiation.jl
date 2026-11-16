# The radiation between the cone, the sapphire window and the sample, from
# first principles (included by LemtaArgon.jl).
#
# - The cone is a gray coil at T_c on the inside of the ISO 5660 frustum
#   (base radius 80 mm, top 40 mm, height 65 mm), 45 mm from the sample. T_c is
#   set so that the cone's radiation transmitted by the window at its starting
#   temperature equals the nominal flux at the sample's centre: the cone was
#   slid in front of the chamber at the start, and its set-point raised to give
#   the nominal flux behind the window (Terrei 2020 thesis, p.54).
# - The window is 5 mm of sapphire (ordinary ray). Its index comes from the
#   one-phonon oscillators and electronic terms of Thomas & Joseph, Johns
#   Hopkins APL Tech. Dig. 9 (1988) 328, Table 1, with the transverse-optic
#   frequencies shifted with temperature as their Fig. 12. Above 1000 cm⁻¹
#   its absorption is the multiphonon edge of their Figs. 6 and 10, read from
#   the scanned figures (±20 %): 0.85 cm⁻¹ at 2000 cm⁻¹ and 0.045 at 2500 cm⁻¹
#   at room temperature, an exponential in wavenumber, growing with temperature
#   as their Fig. 10. A plane slab: Fresnel reflection at both faces, internal
#   attenuation, absorptance = 1 − transmittance − reflectance.
# - The window's temperature follows its own balance at the sample's axis
#   (no conduction along it): absorption of the cone's and the sample's
#   radiation, emission from both faces, natural convection on both faces.
# - The sample's face exchanges with the window as two parallel planes: it
#   receives the cone's transmitted radiation, the window's emission and its
#   own emission reflected back by the window (lagged by one iteration).
#
# All spectral sums run over wavenumber ν [cm⁻¹].

const C1_PLANCK = 3.741771852e-8   # W/(m² (cm⁻¹)⁴): 2πhc², ν in cm⁻¹
const C2_PLANCK = 1.438776877      # cm K
const σ_PLANCK = 5.670374419e-8
const NU = collect(20.0:4.0:25000.0)   # cm⁻¹ (0.4–500 µm)
const DNU = 4.0

planck(ν, T) = C1_PLANCK * ν^3 / expm1(C2_PLANCK * ν / T)

# ∑ f(ν) E_b(ν, T) dν over NU for a spectral vector f
band(f::AbstractVector, T) = DNU * sum(f[i] * planck(NU[i], T) for i in eachindex(NU))

# Sapphire, ordinary ray (Thomas & Joseph 1988, Table 1): ν_TO [cm⁻¹], Δε, γ/ν_TO,
# and the temperature slope of ν_TO [cm⁻¹/K] (their Fig. 12; none given for 806).
const SAPPHIRE_MODES = (
    (385.0, 0.30, 0.015, 1.26e-2),
    (442.0, 2.70, 0.010, 1.06e-2),
    (569.0, 3.00, 0.020, 1.35e-2),
    (635.0, 0.30, 0.020, 1.46e-2),
    (806.0, 0.004, 0.035, 0.0),
)
const SAPPHIRE_EPS_INF = 1 + 1.058 + 1.024   # the two electronic terms, far below them

# Multiphonon absorption coefficient [cm⁻¹] at ν ≥ 1000 cm⁻¹ and T [K].
function sapphire_multiphonon(ν, T)
    a = 0.0014 + 0.0014 * (ν - 2000) / 1000   # 1/K, Fig. 10 at 2000 and 2500 cm⁻¹
    return 0.85 * exp(-(ν - 2000) / 161) * exp(max(a, 0.0) * (T - 300))
end

# Complex index (n, k) of sapphire at ν [cm⁻¹] and T [K].
function sapphire_nk(ν, T)
    ε = complex(SAPPHIRE_EPS_INF)
    for (ν0, Δε, g, s) in SAPPHIRE_MODES
        νj = ν0 - s * (T - 295)
        ε += Δε * νj^2 / (νj^2 - ν^2 - im * ν * g * νj)
    end
    N = sqrt(ε)
    n, k = real(N), imag(N)
    ν >= 1000 && (k = sapphire_multiphonon(ν, T) / (4π * ν))
    return n, k
end

"""
Normal transmittance, reflectance and absorptance of a sapphire plate of
thickness `d` [m] at temperature `T` [K], as spectra over `NU`.
"""
function sapphire_plate(T; d = 0.005)
    τ = similar(NU)
    ρ = similar(NU)
    α = similar(NU)
    for (i, ν) in enumerate(NU)
        n, k = sapphire_nk(ν, T)
        R = ((n - 1)^2 + k^2) / ((n + 1)^2 + k^2)
        t = exp(-4π * ν * k * 100d)
        den = 1 - R^2 * t^2
        τ[i] = (1 - R)^2 * t / den
        ρ[i] = R + R * (1 - R)^2 * t^2 / den
        α[i] = max(1 - τ[i] - ρ[i], 0.0)
    end
    return (; τ, ρ, α)
end

# The plate tabulated every 25 K, interpolated linearly in temperature.
const PLATE_T = collect(250.0:25.0:1100.0)
const PLATE = Ref{Any}(nothing)
function plate(T)
    PLATE[] === nothing && (PLATE[] = [sapphire_plate(t) for t in PLATE_T])
    tab = PLATE[]
    Tc = clamp(T, PLATE_T[1], PLATE_T[end])
    i = min(searchsortedlast(PLATE_T, Tc), length(PLATE_T) - 1)
    w = (Tc - PLATE_T[i]) / (PLATE_T[i+1] - PLATE_T[i])
    return (
        τ = (1 - w) .* tab[i].τ .+ w .* tab[i+1].τ,
        ρ = (1 - w) .* tab[i].ρ .+ w .* tab[i+1].ρ,
        α = (1 - w) .* tab[i].α .+ w .* tab[i+1].α,
    )
end

"""
View factor from a point on the axis, facing the cone, at distance `z` [m]
below the base of the ISO 5660 frustum, to its inside wall: what enters the
base opening minus what leaves by the top opening.
"""
cone_view(z; R_base = 0.080, R_top = 0.040, H = 0.065) =
    R_base^2 / (R_base^2 + z^2) - R_top^2 / (R_top^2 + (z + H)^2)

# Sapphire heat capacity [J/(kg K)]: corundum, NIST Shomate (298–2300 K).
function c_sapphire(T)
    t = T / 1000
    return (102.4290 + 38.74980t - 15.91090t^2 + 2.628181t^3 - 3.007551 / t^2) / 0.1019613
end

"""
    ChamberGeometry(; cone_gap = 0.045, window_gap = 0.020, thickness = 0.005,
                    ε_cone = 0.9, h_window = 7.0, T_room = 293.15,
                    F_window_sample = 0.9)

The cone to sample distance, the gap between the sample and the window's inner
face, the window's thickness, the coil's emissivity, the natural convection on
each face of the window, the room temperature on the cone's side, and the
share of the window's view at its centre taken by the sample and its holder.
"""
Base.@kwdef struct ChamberGeometry
    cone_gap::Float64 = 0.045
    window_gap::Float64 = 0.020
    thickness::Float64 = 0.005
    ε_cone::Float64 = 0.9
    h_window::Float64 = 7.0
    T_room::Float64 = 293.15
    F_window_sample::Float64 = 0.9
end

"""
The coil temperature [K] that sends `q` [W/m²] through the window at `T_w0`
to the sample's centre.
"""
function cone_temperature(q, g::ChamberGeometry; T_w0 = g.T_room)
    F = cone_view(g.cone_gap)
    τ = plate(T_w0).τ
    lo, hi = 500.0, 2000.0
    for _ = 1:60
        T = (lo + hi) / 2
        F * g.ε_cone * band(τ, T) < q ? (lo = T) : (hi = T)
    end
    return (lo + hi) / 2
end

"""
    window_history(q, T_s, T_a, t_end, g; T_w0, dt = 1.0) -> NamedTuple

The window's temperature and what it sends to the sample, at each second to
`t_end`, for nominal flux `q` [W/m²], the sample's surface temperature `T_s(t)`
[K] and emissivity `ε_s(T)`, and the chamber's gas temperature `T_a(t)` [K]:
`T_w`, the cone's radiation transmitted to the sample `q_cone`, and the
long-wave irradiation of the sample `G` (the window's emission, the sample's
own emission reflected back, the chamber walls at T_a beyond the window) [W/m²].
"""
function window_history(
    q,
    T_s,
    ε_s,
    T_a,
    t_end,
    g::ChamberGeometry;
    T_w0 = g.T_room,
    dt = 1.0,
)
    T_c = cone_temperature(q, g; T_w0)
    F_s = cone_view(g.cone_gap)
    F_w = cone_view(g.cone_gap - g.window_gap - g.thickness)   # window's front face
    E_c = [g.ε_cone * planck(ν, T_c) for ν in NU]
    ρc = 3980.0 * g.thickness
    ts = collect(0.0:dt:Float64(t_end))
    Tw = similar(ts)
    qc = similar(ts)
    G = similar(ts)
    Tw[1] = T_w0
    F_ws = g.F_window_sample
    for (k, t) in enumerate(ts)
        p = plate(Tw[k])
        Ts, es = T_s(t), ε_s(T_s(t))
        Ta = T_a(t)
        # sample side: cone through the window, window emission, reflection of
        # the sample's own emission
        qc[k] = F_s * DNU * sum(p.τ[i] * E_c[i] for i in eachindex(NU))
        G[k] = band(p.α, Tw[k]) + es * band(p.ρ, Ts)
        k == length(ts) && break
        # window: absorbed − emitted − convected, per unit area
        absorbed =
            F_w * DNU * sum(p.α[i] * E_c[i] for i in eachindex(NU)) +
            (1 - F_w) * band(p.α, g.T_room) +
            F_ws * es * band(p.α, Ts) +
            (1 - F_ws) * band(p.α, Ta)
        emitted = 2 * band(p.α, Tw[k])
        convected = g.h_window * (Tw[k] - g.T_room) + g.h_window * (Tw[k] - Ta)
        Tw[k+1] = Tw[k] + dt * (absorbed - emitted - convected) / (ρc * c_sapphire(Tw[k]))
    end
    return (; t = ts, T_w = Tw, q_cone = qc, G, T_c)
end

"""
Absorptivity of a surface with spectral absorptivity `a(ν)` for the cone's
radiation at `T_c` transmitted by the window at `T_w`.
"""
function absorptivity_behind_window(a::AbstractVector, T_c, T_w)
    τ = plate(T_w).τ
    w = [τ[i] * planck(NU[i], T_c) for i in eachindex(NU)]
    return sum(a .* w) / sum(w)
end

"""
The window chamber of the decks (chamber.toml [window]) as a `ChamberGeometry`.
"""
function chamber_geometry(ch)
    w = ch["window"]
    return ChamberGeometry(
        cone_gap = w["cone_gap"],
        window_gap = w["window_gap"],
        thickness = w["thickness"],
        ε_cone = w["cone_emissivity"],
        h_window = w["h_window"],
        F_window_sample = w["view_window_sample"],
    )
end

"""
Spectral absorptivity over `NU` of the calcium silicate: two bands fitted to
its Planck-mean emissivity (chamber.toml [silcal]).
"""
function silcal_spectrum(ch)
    p = ch["silcal"]["bands"]
    return [ν < p["edge"] ? p["long"] : p["short"] for ν in NU]
end

"""
Spectral absorptivity over `NU` of the virgin wood of `species`
(decks/wood_absorptivity.csv; fir takes spruce's), held at its end values
outside 0.5–22.25 µm.
"""
function wood_spectrum(species)
    lines = filter(
        l -> !startswith(l, "#"),
        readlines(joinpath(DECKS, "wood_absorptivity.csv")),
    )
    header = split(lines[1], ',')
    j = findfirst(==(species == "fir" ? "spruce" : species), header)
    λ = [parse(Float64, split(l, ',')[1]) for l in lines[2:end]]
    a = [parse(Float64, split(l, ',')[j]) for l in lines[2:end]]
    return [interp(λ, a, 1e4 / ν) for ν in NU]
end

"""
    window_series(build, q, T_a, ε_s, t_end, g; passes = 3, reltol, abstol, dtmax)
        -> series

The window's history under a sample: `build(series)` returns the sample's
`(problem, nodes)` for a window history `series` (`window_history`). The
window is computed from the sample's surface temperature of the previous solve
(the first from a surface warming to T_a + 400 K over a minute), `passes`
times; the last history is returned, with the return codes of the passes'
solves (`retcodes`).
"""
function window_series(
    build,
    q,
    T_a,
    ε_s,
    t_end,
    g;
    passes::Int = 3,
    reltol = 1e-6,
    abstol = 1e-6,
    dtmax = 5.0,
    use_ale::Bool = false,
)
    T_s = t -> T_a(t) + 400 * (1 - exp(-t / 60))
    series = window_history(q, T_s, ε_s, T_a, t_end, g)
    retcodes = String[]
    for _ = 1:passes
        prob, _ = build(series)
        sol = solve(prob; reltol, abstol, dtmax, use_ale)
        push!(retcodes, string(sol.retcode))
        ts, Ts = copy(sol.t), copy(sol.surface_T)
        T_s = t -> interp(ts, Ts, t)
        series = window_history(q, T_s, ε_s, T_a, t_end, g)
    end
    return merge(series, (; retcodes))
end

# The exposed face under the cone and window: the cone's radiation transmitted
# by the window (absorbed at the surface's absorptivity), the window's
# long-wave irradiation (absorbed and emitted at the surface's emissivity) and
# convection at `h` to T_a.
window_face(series, h, T_a) =
    RadiativeFluxBC(t -> interp(series.t, series.q_cone, t)) +
    RadiativeBC(T_inf = t -> (interp(series.t, series.G, t) / σ_PLANCK)^0.25) +
    ConvectiveBC(h = Float64(h), T_inf = T_a)
