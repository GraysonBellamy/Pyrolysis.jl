# The LEMTA moist-wood case (doi:10.57745/CBLZ8P; Diarra et al., Fire Saf. J.
# 165 (2026) 104951): wood blocks at 0–30 % moisture under a cone heater at
# 4 kW/m² for 2.5 h, then cooling, with temperatures at 12 depths and mass.
# Data loading and the defect list, species decks to Pyrolysis.jl problems,
# the evaporation laws compared, model outputs at the thermocouple depths, and
# the comparison metrics. See README.md.
module LemtaMoistWood

using Pyrolysis
using Printf
using SHA
using Statistics
using TOML

export lemta_data,
    check_data,
    Run,
    load_run,
    load_runs,
    run_names,
    deck,
    initial_water_g,
    material,
    problem,
    simulate,
    write_series,
    read_series,
    series_path,
    RESULTS,
    LAWS,
    VARIANTS,
    k_across,
    nrms,
    nrms_by_depth,
    mass_loss,
    crossing_minutes,
    rate_ratio,
    interp,
    p_sat,
    heat_of_desorption

const HERE = @__DIR__
const DECKS = joinpath(HERE, "decks")
const RESULTS = joinpath(HERE, "results")
const CALIBRATION = joinpath(HERE, "calibration.toml")

const AREA = 0.01               # m², the 10 cm × 10 cm exposed face
const THICKNESS = 0.05          # m
const VOLUME = AREA * THICKNESS # m³
const FLUX = 4.0e3              # W/m², incident, set with the fluxmeter at the sample centre
const P_ATM = 101325.0          # Pa
const SWITCH = 3.0              # s, width of the smooth step when the cone is withdrawn

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

lemta_data() =
    get(ENV, "LEMTA_DATA", joinpath(homedir(), "Documents", "data", "lemta-CBLZ8P"))

"""Check every file of the dataset against its pinned SHA-256 (data_files.toml)."""
function check_data(; root = lemta_data())
    isdir(root) ||
        error("LEMTA data not found at $root (set LEMTA_DATA, or run fetch_data.sh)")
    pinned = TOML.parsefile(joinpath(HERE, "data_files.toml"))["files"]
    for f in pinned
        path = joinpath(root, f["name"])
        isfile(path) || error("missing $(f["name"]) in $root (run fetch_data.sh)")
        got = bytes2hex(open(sha256, path))
        got == f["sha256"] ||
            error("$(f["name"]) differs from the published file (SHA-256 $got)")
    end
    return root
end

# Runs the published data cannot support as recorded.
#   :excluded   oak_test2_0%: the channels are not ordered with depth (12 mm at
#               91 °C, 14 mm at 186 °C, 30 mm at 151 °C at 2.5 h); scrambled.
#   :loss_only  the two fir files that start at 535–554 g against a dry mass of
#               about 217 g: the holder was weighed with the sample, so only the
#               mass lost is read, and the initial water is MC × ρ_dry × V.
const DEFECTS = Dict(
    "oak_test2_0%_511" => :excluded,
    "fir_17,67%_433" => :loss_only,
    "fir_9,1%_433" => :loss_only,
)

"""
One test as recorded. `T` is in K, one column per thermocouple at `depths` [m]
from the heated face; `mass` [g] is empty for the dry runs. `t_off` is when the
cone was withdrawn (Inf when the record ends under the cone). `T_amb` is the
shallowest thermocouple's first reading, and `T0` the first row, which is the
initial profile of every run (`spruce_test2_0%` starts with a 20 → 40 °C
gradient).
"""
struct Run
    name::String
    species::String
    mc::Float64
    rho::Float64
    repeat::Int
    t::Vector{Float64}
    depths::Vector{Float64}
    T::Matrix{Float64}
    mass::Vector{Float64}
    t_off::Float64
    T_amb::Float64
    T0::Vector{Float64}
    defect::Symbol
end

Base.show(io::IO, r::Run) = print(
    io,
    "Run($(r.name): $(r.species), MC $(round(100r.mc, digits = 2)) %, ρ $(r.rho), ",
    "$(length(r.depths)) TCs, $(round(r.t[end] / 3600, digits = 2)) h)",
)

function parse_name(name)
    m = match(r"^([a-z]+)_(?:test(\d)_)?([\d,]+)%_(\d+)$", name)
    m === nothing && error("unexpected file name $name")
    species = m.captures[1]
    rep = m.captures[2] === nothing ? 0 : parse(Int, m.captures[2])
    mc = parse(Float64, replace(m.captures[3], "," => ".")) / 100
    rho = parse(Float64, m.captures[4])
    return species, mc, rho, rep
end

run_names(; root = lemta_data()) =
    sort([splitext(f)[1] for f in readdir(root) if endswith(f, ".dat")])

"""
    load_run(name; root) -> Run

Parse one file by its own header: the thermocouple columns present (`2mm` …
`40mm`, some missing), the mass in g (`PM (g)`, written `PM (s)` in two
files; absent in the dry runs), LF or CRLF endings, `Temps` or `Time`.
"""
function load_run(name::AbstractString; root = lemta_data())
    species, mc, rho, rep = parse_name(name)
    lines = readlines(joinpath(root, name * ".dat"))
    header = strip.(split(chomp(replace(lines[1], '\r' => "")), '\t'))
    itime = findfirst(h -> occursin(r"^(Temps|Time) \(s\)$", h), header)
    itime === nothing && error("$name: no time column in $header")
    tcs = [
        (i, parse(Int, h[1:(end-2)])) for
        (i, h) in enumerate(header) if occursin(r"^\d+mm$", h)
    ]
    imass = findfirst(h -> h in ("PM (g)", "PM (s)"), header)
    rows = Vector{Vector{Float64}}()
    for l in lines[2:end]
        f = split(strip(l), '\t')
        length(f) < length(header) && continue
        push!(rows, parse.(Float64, f))
    end
    t = [r[itime] for r in rows]
    T = [r[i] + 273.15 for r in rows, (i, _) in tcs]
    mass = imass === nothing ? Float64[] : [r[imass] for r in rows]
    depths = [d / 1000 for (_, d) in tcs]
    defect = get(DEFECTS, name, :none)
    T_amb = mean(T[1:min(5, end), 1])
    return Run(
        name,
        species,
        mc,
        rho,
        rep,
        t,
        depths,
        T,
        mass,
        cone_off(t, T[:, 1]),
        T_amb,
        T[1, :],
        defect,
    )
end

# The cone's withdrawal: the first time after 2.3 h at which the shallowest
# thermocouple falls by more than 6 K within the next 30 s (the withdrawal
# drops it by 10–37 K in 30 s; noise spikes reach 2–3 K), less the 2 mm
# thermocouple's lag, taken as the first second of that fall.
function cone_off(t, Ts)
    for i in eachindex(t)
        t[i] < 8300 && continue
        j = searchsortedfirst(t, t[i] + 30)
        j > length(t) && break
        if Ts[j] - Ts[i] < -6
            k = i
            while k < j && Ts[k+1] - Ts[k] > -0.5
                k += 1
            end
            return t[k]
        end
    end
    return Inf
end

"""All runs, the excluded one dropped unless `excluded = true`."""
function load_runs(; root = lemta_data(), excluded = false)
    runs = [load_run(n; root) for n in run_names(; root)]
    excluded || filter!(r -> r.defect != :excluded, runs)
    return runs
end

"""
Initial water [g]: MC/(1 + MC) × m₀ (dry basis, Diarra et al. Eq. 1) from the
run's own first mass, or MC ρ_dry V when the first mass is not the sample's.
"""
function initial_water_g(r::Run)
    r.mc == 0 && return 0.0
    r.defect == :loss_only && return r.mc * r.rho * VOLUME * 1000
    return r.mc / (1 + r.mc) * r.mass[1]
end

# ---------------------------------------------------------------------------
# Decks and properties
# ---------------------------------------------------------------------------

deck(species::AbstractString) = TOML.parsefile(joinpath(DECKS, species * ".toml"))

# Piecewise linear, held constant outside the table.
function interp(xs::AbstractVector, ys::AbstractVector, x::Real)
    x <= xs[1] && return float(ys[1])
    x >= xs[end] && return float(ys[end])
    i = searchsortedlast(xs, x)
    w = (x - xs[i]) / (xs[i+1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i+1]
end

"""
Conductivity across the grain [W/(m K)], Flity et al., Int. J. Therm. Sci. 195
(2024) 108661, Eq. 9: one law of dry density and temperature for every
species, fitted to hot-wire data at 20–160 °C (rms 0.0065 W/(m K)).
"""
k_across(ρ_dry, T) = 0.0386 + 0.1708 * (ρ_dry / 1000) * sqrt(T / 433)

"""
Heat capacity of dry wood [J/(kg K)], Flity et al. 2025 §4.2 (the mean of
Flity et al. 2024 Table 2 for eight species, 20–90 °C, extrapolated above).
"""
c_wood(T) = 1131 + 4.67 * (T - 273)

# Liquid water, Shomate (Chase 1998), as in examples/07.
function c_water(T)
    t = T / 1000
    return (-203.6060 + 1523.290t - 3196.413t^2 + 2474.455t^3 + 3.855326 / t^2) / 0.01801528
end
c_vapour(T) = 1716.7 + 0.4154T + 1.530e-4T^2   # examples/07, NIST fit 300–1200 K

"""Saturation pressure of water [Pa] (Antoine, 1–374 °C)."""
function p_sat(T)
    θ = T - 273.15
    a, b, c = θ < 100 ? (8.07131, 1730.63, 233.426) : (8.14019, 1810.94, 244.485)
    return 133.322 * 10^(a - b / (c + θ))
end

"""
Heat taken by desorbing one kg of bound water at `T` [K] with `X` kg of water
per kg of dry wood [J/kg]: the latent heat of free water at `T` (linear,
2501 kJ/kg at 0 °C, 2257 at 100 °C) plus the differential heat of sorption of
examples/07, 875.6 exp(−11.22 X) kJ/kg, scaled by `sorption`. Heats of
reaction are quoted at the local temperature (TRG §1.5.3), hence L(T).
"""
heat_of_desorption(T, X; sorption = 1.0) =
    2.501e6 - 2369.0 * (T - 273.15) +
    sorption * 875.6441221735336e3 * exp(-11.21684273316275X)

"""
The evaporation laws compared, as (A [1/s], E [J/mol], T_min [K], gate width [K]):
`:arrhenius` is examples/07's set; `:gated` the same switched on at 373.15 K
over 2 K (kinetics audit §4.1a); `:chan` the drying set the literature quotes
from Chan, Kelbon & Krieger 1985 (kinetics audit §4.1). `:onset343` …
`:onset363` and `:onset383` switch examples/07's set on at 70, 80, 90 and 110 °C
over 5 K, and `:onset373w10` at 100 °C over 10 K: a scan of where the latent
sink sits, not candidate laws.
"""
const LAWS = Dict(
    :arrhenius => (A = 2.997e5, E = 5.0e4, T_min = 0.0, width = 1.0),
    :gated => (A = 2.997e5, E = 5.0e4, T_min = 373.15, width = 2.0),
    :chan => (A = 5.13e10, E = 88.0e3, T_min = 0.0, width = 1.0),
    # exploratory: examples/07's set switched on below the boiling point
    :onset343 => (A = 2.997e5, E = 5.0e4, T_min = 343.15, width = 5.0),
    :onset353 => (A = 2.997e5, E = 5.0e4, T_min = 353.15, width = 5.0),
    :onset363 => (A = 2.997e5, E = 5.0e4, T_min = 363.15, width = 5.0),
    :onset373w10 => (A = 2.997e5, E = 5.0e4, T_min = 373.15, width = 10.0),
    :onset383 => (A = 2.997e5, E = 5.0e4, T_min = 383.15, width = 5.0),
)

const IWOOD, IWATER, IVAPOUR, IAIR = 1, 2, 3, 4

"""
    material(d; law, κ, λ = 5e-7, sorption = 1.0, k_rule = PARALLEL, A_scale = 1.0,
             k_scale = 1.0)

The deck's wood. `law = :dry` gives the dry wood alone (thermal model only);
otherwise dry wood, bound water, water vapour and pore air, with the water
desorbing by the evaporation law `law` (a key of `LAWS`). The dry wood's
volume share is one at its dry density; the water holds its own volume at
1000 kg/m³, and under `k_rule = PARALLEL` adds to the conductivity by that
share at 0.6 W/(m K), which gives the +8 to +26 % that Flity et al. 2024
(Table 6) measured at about 10 % moisture across the grain; `SERIES` removes
the effect. `κ` [m²] and `λ` [m²/s] are the dry wood's Darcy permeability and
lumped vapour transfer coefficient; the water blocks neither and carries no
gas. `A_scale` multiplies the law's pre-exponential factor and `k_scale` the
dry wood's conductivity.
"""
function material(
    d;
    law::Symbol,
    κ::Real = 1e-15,
    λ::Real = 5e-7,
    sorption::Real = 1.0,
    k_rule = PARALLEL,
    A_scale::Real = 1.0,
    k_scale::Real = 1.0,
)
    ρ = Float64(d["rho_dry"])
    Tε = Float64.(d["emissivity_T_C"]) .+ 273.15
    ε_tab = Float64.(d["emissivity"])
    εlaw = CallableProperty(T -> interp(Tε, ε_tab, T))
    wood = SolidComponent(
        :WOOD;
        ρ = ρ,
        c = CallableProperty(c_wood),
        k = CallableProperty(T -> k_scale * k_across(ρ, T)),
        ε = εlaw,
        λ = λ,
        κ = κ,
        ρ_intrinsic = 1500.0,       # cell-wall substance (Siau)
    )
    law == :dry &&
        return Material(name = Symbol(d["species"]), components = (wood,), reactions = ())
    water = LiquidComponent(
        :MOISTURE;
        ρ = 1000.0,
        c = CallableProperty(c_water),
        k = 0.6,
        ε = εlaw,
    )
    vapour = GaseousComponent(
        :WATER_VAPOR;
        M = 0.018015,
        c = CallableProperty(c_vapour),
        k = 0.025,
    )
    air = GaseousComponent(:AIR; M = 0.02897, c = 1006.0, k = 0.026)
    p = LAWS[law]
    h = StateDependentProperty(
        (T, ξ) -> heat_of_desorption(T, ξ[IWATER] / max(ξ[IWOOD], 1e-12); sorption),
    )
    evaporation = Reaction(
        :evaporation,
        IWATER => IVAPOUR;
        A = p.A * A_scale,
        E = p.E,
        h = h,
        n = 1.0,
        T_min = p.T_min,
        gate_width = p.width,
    )
    return Material(
        name = Symbol(d["species"]),
        components = (wood, water, vapour, air),
        reactions = (evaporation,),
        mixing = MaterialMixing(
            k = MixingSpec(k_rule),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
    )
end

# ---------------------------------------------------------------------------
# Scenario
# ---------------------------------------------------------------------------

# 1 under the cone, 0 after its withdrawal at `t_off`, a smooth step of width SWITCH.
cone_on(t, t_off) = isfinite(t_off) ? 0.5 * (1 - tanh((t - t_off) / SWITCH)) : 1.0

"""
    problem(r::Run; law, cells = 100, t_end, flux = FLUX, h_scale = 1.0,
            rear_radiation = false, absorptivity, h_front, h_rear, κ, κ_scale,
            λ, sorption, k_rule, A_scale, k_scale) -> (problem, mesh nodes)

The 1D slab of 50 mm across the grain: the exposed face (`:top`, z = 50 mm)
absorbs `absorptivity × flux` while the cone is on and re-radiates at ε(T_s)
to the ambient, with convection at the deck's heating and cooling
coefficients; the back (`:bottom`) loses heat at the deck's `h_rear` (the
whole loss Flity et al. 2025 fitted there) plus radiation when
`rear_radiation`. Both faces are open to the ambient air. The initial
temperature is the run's first row, interpolated in depth. `κ` defaults to
the deck's `permeability`, scaled by `κ_scale`; the other keywords are those
of [`material`](@ref) or override the deck. With `calibrated`, the incident
flux and the conductivity take the species' factors fitted to its dry runs
(calibration.toml; a species without dry runs keeps the deck).
"""
function problem(
    r::Run;
    law::Symbol,
    cells::Int = 100,
    t_end::Real = r.t[end],
    flux::Real = FLUX,
    absorptivity::Union{Nothing,Real} = nothing,
    h_front::Union{Nothing,Tuple} = nothing,
    h_rear::Union{Nothing,Real} = nothing,
    rear_radiation::Bool = false,
    κ::Union{Nothing,Real} = nothing,
    κ_scale::Real = 1.0,
    λ::Real = 5e-7,
    sorption::Real = 1.0,
    k_rule = PARALLEL,
    A_scale::Real = 1.0,
    k_scale::Real = 1.0,
    calibrated::Bool = false,
    save_every::Real = 10.0,
)
    d = deck(r.species)
    if calibrated
        c = get(
            TOML.parsefile(CALIBRATION),
            r.species,
            Dict("flux_factor" => 1.0, "k_factor" => 1.0),
        )
        flux *= c["flux_factor"]
        k_scale *= c["k_factor"]
    end
    κ = κ_scale * something(κ, Float64(get(d, "permeability", 1e-15)))
    mat =
        material(d; law = r.mc == 0 ? :dry : law, κ, λ, sorption, k_rule, A_scale, k_scale)
    α = something(absorptivity, Float64(d["absorptivity"]))
    hh, hc =
        something(h_front, (Float64(d["h_front_heating"]), Float64(d["h_front_cooling"])))
    hr = something(h_rear, Float64(d["h_rear"]))
    Tε = Float64.(d["emissivity_T_C"]) .+ 273.15
    ε_tab = Float64.(d["emissivity"])
    εs(T) = interp(Tε, ε_tab, T)
    Ta = r.T_amb
    t_off = r.t_off

    top =
        RadiativeFluxBC(t -> flux * cone_on(t, t_off); absorptivity = α) +
        RadiativeBC(ε = εs, T_inf = Ta) +
        ConvectiveBC(h = t -> hc + (hh - hc) * cone_on(t, t_off), T_inf = Ta)
    bottom = ConvectiveBC(h = hr, T_inf = Ta)
    rear_radiation && (bottom = bottom + RadiativeBC(ε = εs, T_inf = Ta))

    nodes = collect(range(0.0, THICKNESS, length = cells + 1))
    mesh = line_mesh(nodes)
    order = sortperm(r.depths)
    dep, T0 = r.depths[order], r.T0[order]
    Tinit = z -> interp(dep, T0, THICKNESS - z)

    if r.mc == 0
        bcs = BoundaryConditionSet(thermal = (top = top, bottom = bottom))
        ξinit = [z -> Float64(d["rho_dry"])]
    else
        reservoir = GasReservoir(P = P_ATM, T = Ta, mole_fractions = (AIR = 1.0,))
        bcs = BoundaryConditionSet(
            thermal = (top = top, bottom = bottom),
            mass = (top = AmbientGasBC(reservoir), bottom = AmbientGasBC(reservoir)),
        )
        ρw = initial_water_g(r) / 1000 / VOLUME
        ξ0 = (Float64(d["rho_dry"]), ρw, 0.0, 0.0)
        ξinit = [
            z -> ξ0[1],
            z -> ξ0[2],
            z -> 0.0,
            z -> background_gas_concentration(mat, :AIR, ξ0, Tinit(z)),
        ]
    end
    prob = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bcs,
        T_initial = Tinit,
        ξ_initial = ξinit,
        tspan = (0.0, Float64(t_end)),
        dt_initial = 1e-2,
        output_times = collect(0.0:save_every:Float64(t_end)),
    )
    return prob, nodes
end

"""
    simulate(r::Run; law, kwargs...) -> NamedTuple

Solve the run's problem and return the series the comparison reads: the
thermocouple temperatures [K] (cell centres interpolated linearly to each
depth, the depths scaled by `probe_scale`), the sample mass [g] (dry wood and
water), the water left [g], the
largest pore pressure and vapour-to-saturation ratio, and the evaporation
rate [g/s] with the share of it taking place where the pore vapour pressure
exceeds the saturation pressure at the cell temperature (where a
vapour-pressure-limited law would not evaporate).
"""
function simulate(
    r::Run;
    law::Symbol = :arrhenius,
    reltol = 1e-6,
    dtmax = 30.0,
    probe_scale = 1.0,
    kwargs...,
)
    prob, nodes = problem(r; law, kwargs...)
    sol = solve(prob; reltol, abstol = 1e-6, dtmax, radiation_model = SURFACE_ABSORPTION)
    mat = prob.material
    nt = length(sol.t)
    centres = 0.5 .* (nodes[1:(end-1)] .+ nodes[2:end])
    widths = diff(nodes)
    zs = THICKNESS .- probe_scale .* r.depths
    Ttc = Matrix{Float64}(undef, nt, length(zs))
    for k = 1:nt, (j, z) in enumerate(zs)
        Ttc[k, j] = interp(centres, view(sol.T, :, k), z)
    end
    NC = size(sol.ξ, 1)
    moist = NC > 1
    mass = zeros(nt)
    water = zeros(nt)
    pmax = fill(P_ATM, nt)
    smax = zeros(nt)
    evap = zeros(nt)
    evap_super = zeros(nt)
    rates = zeros(max(length(mat.reactions), 1))
    heats = similar(rates)
    for k = 1:nt
        for c in eachindex(widths)
            ξ = ntuple(j -> sol.ξ[j, c, k], NC)
            mass[k] += ξ[1] * widths[c] * AREA * 1000
            moist || continue
            mass[k] += ξ[IWATER] * widths[c] * AREA * 1000
            water[k] += ξ[IWATER] * widths[c] * AREA * 1000
            T = sol.T[c, k]
            pore = Pyrolysis.Physics.pore_state(ξ, T, mat)
            pmax[k] = max(pmax[k], pore.P)
            s = pore.p[IVAPOUR] / p_sat(T)
            smax[k] = max(smax[k], s)
            Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, mat)
            e = rates[1] * widths[c] * AREA * 1000
            evap[k] += e
            s > 1 && (evap_super[k] += e)
        end
    end
    return (;
        t = sol.t,
        depths = r.depths,
        T = Ttc,
        T_surface = sol.surface_T,
        T_back = sol.extras.back_temperature,
        mass,
        water,
        p_max = pmax,
        supersaturation = smax,
        evaporation = evap,
        evaporation_supersaturated = evap_super,
        retcode = string(sol.retcode),
        solve_time = sol.solve_time,
    )
end

# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------

# Conductivity across the grain at 20 °C that Diarra et al. Table 1 quotes for
# each species (Flity et al. 2024's samples, at their own densities).
const K_TABLE1 = Dict(
    "larch" => 0.123,
    "oak" => 0.125,
    "fir" => 0.095,
    "spruce" => 0.096,
    "poplar" => 0.072,
)

"""
The variants run, by name: `(runs = :dry | :moist, kwargs = r -> NamedTuple)`,
the keywords handed to [`simulate`](@ref) for run `r`.
"""
const VARIANTS = Dict{String,NamedTuple}(
    # dry runs: the thermal model and the diagnostics of its residual
    "dry" => (runs = :dry, kwargs = r -> (;)),
    "dry_cells200" => (runs = :dry, kwargs = r -> (cells = 200,)),
    "dry_tight" => (runs = :dry, kwargs = r -> (reltol = 1e-8, dtmax = 5.0)),
    "dry_handoff_bc" => (
        runs = :dry,
        kwargs = r -> (h_front = (9.5, 9.5), h_rear = 3.5, rear_radiation = true),
    ),
    "dry_probe_flity" => (runs = :dry, kwargs = r -> (probe_scale = 1.035,)),
    "dry_k_table1" => (
        runs = :dry,
        kwargs = r -> (k_scale = K_TABLE1[r.species] / k_across(r.rho, 293.15),),
    ),
    "dry_flux3.8" => (runs = :dry, kwargs = r -> (flux = 3.8e3,)),
)
VARIANTS["dry_cal"] = (runs = :dry, kwargs = r -> (calibrated = true,))
# Moist runs: the deck as published (`<law>`) and with the dry runs' factors
# (`<law>_cal`); the sensitivities on the calibrated deck.
for law in (:arrhenius, :gated, :chan)
    l = String(law)
    VARIANTS[l] = (runs = :moist, kwargs = r -> (; law))
    VARIANTS[l*"_cal"] = (runs = :moist, kwargs = r -> (; law, calibrated = true))
    for (tag, kw) in (
        "kappa0.1" => (κ_scale = 0.1,),
        "kappa10" => (κ_scale = 10.0,),
        "vent" => (κ = 1e-15,),
        "sorption0" => (sorption = 0.0,),
        "sorption2" => (sorption = 2.0,),
        "kseries" => (k_rule = SERIES,),
        "cells200" => (cells = 200,),
        "tight" => (reltol = 1e-8, dtmax = 10.0),
    )
        VARIANTS[l*"_cal_"*tag] =
            (runs = :moist, kwargs = r -> (; law, calibrated = true, kw...))
    end
end
for law in (:onset343, :onset353, :onset363, :onset373w10, :onset383)
    VARIANTS[String(law)*"_cal"] = (runs = :moist, kwargs = r -> (; law, calibrated = true))
end

# ---------------------------------------------------------------------------
# Series files
# ---------------------------------------------------------------------------

series_path(variant, name) = joinpath(RESULTS, variant, name * ".csv")

function write_series(path, s)
    mkpath(dirname(path))
    open(path, "w") do io
        tcs = [@sprintf("T_%gmm_K", 1000d) for d in s.depths]
        println(
            io,
            join(
                [
                    "time_s";
                    tcs;
                    [
                        "T_surface_K",
                        "T_back_K",
                        "mass_g",
                        "water_g",
                        "p_max_Pa",
                        "pv_over_psat_max",
                        "evaporation_g_s",
                        "evaporation_supersaturated_g_s",
                    ]
                ],
                ",",
            ),
        )
        for k in eachindex(s.t)
            vals = [
                s.t[k];
                s.T[k, :];
                [
                    s.T_surface[k],
                    s.T_back[k],
                    s.mass[k],
                    s.water[k],
                    s.p_max[k],
                    s.supersaturation[k],
                    s.evaporation[k],
                    s.evaporation_supersaturated[k],
                ]
            ]
            println(io, join((@sprintf("%.6g", v) for v in vals), ","))
        end
        println(io, "# retcode=$(s.retcode) solve_time=$(round(s.solve_time, digits = 2))")
    end
    return path
end

function read_series(path)
    lines = filter(l -> !startswith(l, "#"), readlines(path))
    header = split(lines[1], ',')
    data = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
    col(name) = [row[findfirst(==(name), header)] for row in data]
    tcs = [h for h in header if occursin(r"^T_[\d.]+mm_K$", h)]
    depths = [parse(Float64, h[3:(end-4)]) / 1000 for h in tcs]
    T = reduce(hcat, [col(h) for h in tcs])
    meta = filter(l -> startswith(l, "# retcode"), readlines(path))
    return (;
        t = col("time_s"),
        depths,
        T,
        T_surface = col("T_surface_K"),
        T_back = col("T_back_K"),
        mass = col("mass_g"),
        water = col("water_g"),
        p_max = col("p_max_Pa"),
        supersaturation = col("pv_over_psat_max"),
        evaporation = col("evaporation_g_s"),
        evaporation_supersaturated = col("evaporation_supersaturated_g_s"),
        retcode = isempty(meta) ? "" : split(split(meta[1])[2], '=')[2],
    )
end

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

"""
NRMS at one depth, MaCFP style: the rms of model minus measurement over the
measured record (model interpolated to the measured times, both within
`window`), over the measured range max − min.
"""
function nrms(t_exp, y_exp, t_mod, y_mod; window = (0.0, Inf))
    keep = findall(t -> window[1] <= t <= min(window[2], t_mod[end]), t_exp)
    isempty(keep) && return NaN
    y = y_exp[keep]
    m = [interp(t_mod, y_mod, t) for t in t_exp[keep]]
    return sqrt(mean((m .- y) .^ 2)) / (maximum(y) - minimum(y))
end

"""NRMS at every thermocouple depth of the run, as Dict(depth_mm => value)."""
function nrms_by_depth(r::Run, s; window = (0.0, Inf))
    out = Dict{Int,Float64}()
    for (j, d) in enumerate(r.depths)
        js = findfirst(x -> isapprox(x, d; atol = 1e-6), s.depths)
        out[round(Int, 1000d)] = nrms(r.t, view(r.T, :, j), s.t, view(s.T, :, js); window)
    end
    return out
end

"""
Mass lost [g] by time `t`: in the measurement `m₀ − m(t)`, each the median
of the record within 30 s (the balance has single-sample spikes of up to
2 g), or in a model series.
"""
function mass_loss(r::Run, t)
    isempty(r.mass) && return NaN
    window(t0) = median(r.mass[findall(x -> abs(x - t0) <= 30, r.t)])
    return window(0.0) - window(t)
end
mass_loss(s::NamedTuple, t) = s.mass[1] - interp(s.t, s.mass, t)

"""
Minutes to rise from 85 to 115 °C at `depth` [m] within the first `t_max` s
(NaN when 115 °C is not reached), from a (t, T [K]) series.
"""
function crossing_minutes(t, T; lo = 358.15, hi = 388.15, t_max = 9000.0)
    a = findfirst(i -> T[i] >= lo && t[i] <= t_max, eachindex(t))
    b = findfirst(i -> T[i] >= hi && t[i] <= t_max, eachindex(t))
    (a === nothing || b === nothing) && return NaN
    return (t[b] - t[a]) / 60
end

"""
Heating rate near 100 °C over the heating rate near 70 °C, each the
temperature change over the 2 min centred on the first crossing (as
lemta_qc/lemta_dwell.py), from a (t, T [K]) series.
"""
function rate_ratio(t, T; t_max = 9000.0)
    function rate(v)
        k = findfirst(i -> T[i] >= v && t[i] <= t_max, eachindex(t))
        (k === nothing || t[k] < 60 || t[k] + 60 > t[end]) && return NaN
        return (interp(t, T, t[k] + 60) - interp(t, T, t[k] - 60)) / 2
    end
    return rate(373.15) / rate(343.15)
end

end # module
