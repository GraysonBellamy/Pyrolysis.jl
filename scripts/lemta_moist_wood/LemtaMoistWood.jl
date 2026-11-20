# The LEMTA moist-wood case (Diarra et al., Fire Saf. J. 165 (2026) 104951):
# wood blocks at 0–30 % moisture under a cone heater at 4 kW/m² for 2.5 h, then
# cooling, with temperatures at 12 depths and mass, the grain across the heat
# flux (doi:10.57745/CBLZ8P) or along it (the authors' workbooks), and the
# moisture profiles of spruce blocks sliced after 0.5, 1 and 2 h. Data loading
# and the defect list, species decks to Pyrolysis.jl problems, the moisture
# model (the phase change of the bound water and its vapour) and the one-way
# evaporation laws compared with it, model outputs at the thermocouple depths,
# and the comparison metrics. See README.md.
module LemtaMoistWood

using Pyrolysis
using Printf
using SHA
using Statistics
using TOML

export lemta_data,
    lemta_along_data,
    check_data,
    Run,
    load_run,
    load_runs,
    slice_runs,
    run_names,
    deck,
    initial_water_g,
    measured_water_g,
    mass_density,
    material,
    problem,
    simulate,
    write_series,
    read_series,
    series_path,
    write_profiles,
    read_profiles,
    profiles_path,
    RESULTS,
    LAWS,
    PHASE_CHANGE_LAWS,
    ISOTHERMS,
    ONE_WAY,
    VARIANTS,
    k_across,
    nrms,
    nrms_by_depth,
    mass_loss,
    water_fraction_lost,
    crossing_minutes,
    rate_ratio,
    tc_accuracy,
    temperature_errors,
    trapz,
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
const D_VAPOUR_AIR = 2.6e-5     # m²/s, water vapour in air at 25 °C and 1 atm (dusty-gas closure)
const LAMBDA = 1.5e-6           # m²/s, transfer coefficient across the grain (a deck's `lambda` overrides)
const MASS_STEP = 1.07          # s, the balance's sampling interval in the along-grain records
const T_UNRECORDED = 293.15     # K, initial and ambient temperature where no thermocouple recorded them

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

lemta_data() =
    get(ENV, "LEMTA_DATA", joinpath(homedir(), "Documents", "data", "lemta-CBLZ8P"))
lemta_along_data() =
    get(ENV, "LEMTA_ALONG_DATA", joinpath(homedir(), "Documents", "data", "lemta-parallel"))

"""
    check_data(; grain = :across, root)

Check the dataset of the grain direction against its pinned SHA-256: across
the grain every published file (data_files.toml); along it the workbooks as
received and the CSVs convert_along.py writes from them
(data_files_along.toml).
"""
function check_data(;
    grain::Symbol = :across,
    root = grain === :along ? lemta_along_data() : lemta_data(),
)
    if grain === :along
        isdir(root) || error(
            "LEMTA along-grain data not found at $root (set LEMTA_ALONG_DATA; see README.md)",
        )
        manifest = TOML.parsefile(joinpath(HERE, "data_files_along.toml"))
        pinned = [manifest["workbooks"]; manifest["files"]]
        hint = "run convert_along.py"
    else
        isdir(root) ||
            error("LEMTA data not found at $root (set LEMTA_DATA, or run fetch_data.sh)")
        pinned = TOML.parsefile(joinpath(HERE, "data_files.toml"))["files"]
        hint = "run fetch_data.sh"
    end
    for f in pinned
        path = joinpath(root, f["name"])
        isfile(path) || error("missing $(f["name"]) in $root ($hint)")
        got = bytes2hex(open(sha256, path))
        got == f["sha256"] ||
            error("$(f["name"]) differs from the pinned file (SHA-256 $got)")
    end
    return root
end

# Runs the data cannot support as recorded.
#   :excluded   oak_test2_0%: the channels are not ordered with depth (12 mm at
#               91 °C, 14 mm at 186 °C, 30 mm at 151 °C at 2.5 h); scrambled.
#   :loss_only  the two fir files that start at 535–554 g against a dry mass of
#               about 217 g: the holder was weighed with the sample, so only the
#               mass lost is read, and the initial water is MC × ρ_dry × V.
#   :mass_gap   fir_20,2%_para: the balance recorded nothing from 0.93 to 2.50 h;
#               the cone is taken off at 2.5 h, where the record resumes with
#               the loss rate already falling.
const DEFECTS = Dict(
    "oak_test2_0%_511" => :excluded,
    "fir_17,67%_433" => :loss_only,
    "fir_9,1%_433" => :loss_only,
    "fir_20,2%_para" => :mass_gap,
)

"""
One test as recorded. `grain` is `:across` or `:along` the heat flux. `T` is
in K, one column per thermocouple at `depths` [m] from the heated face, at the
times `t` [s]; the along-grain runs at about 20 % have no thermocouples (`T`
has no columns, `t` is the balance's clock). `mass` [g] is empty for the dry
runs and is read at `t_mass` [s], the balance's own clock. `t_off` is when the
cone was withdrawn (Inf when the record ends under the cone). `T_amb` is the
shallowest thermocouple's first reading, and `T0` the first row, which is the
initial profile of every run (`spruce_test2_0%` starts with a 20 → 40 °C
gradient); without thermocouples both are `T_UNRECORDED`.
"""
struct Run
    name::String
    species::String
    grain::Symbol
    mc::Float64
    rho::Float64
    repeat::Int
    t::Vector{Float64}
    depths::Vector{Float64}
    T::Matrix{Float64}
    t_mass::Vector{Float64}
    mass::Vector{Float64}
    t_off::Float64
    T_amb::Float64
    T0::Vector{Float64}
    defect::Symbol
end

Base.show(io::IO, r::Run) = print(
    io,
    "Run($(r.name): $(r.species) $(r.grain), MC $(round(100r.mc, digits = 2)) %, ρ $(r.rho), ",
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

"""
    run_names(; grain = :across, root) -> names

Across the grain the published files, `<species>_[testN_]<MC>%_<ρ>`; along it
the workbooks' sheets, `<species>_<MC>%_para`.
"""
function run_names(;
    grain::Symbol = :across,
    root = grain === :along ? lemta_along_data() : lemta_data(),
)
    grain === :along && return sort(collect(keys(along_files(root))))
    return sort([splitext(f)[1] for f in readdir(root) if endswith(f, ".dat")])
end

# The along-grain runs by name, each with its CSV (one sheet of a workbook,
# `<Species>_[Pp]ara__<MC>%.csv`, the decimal written either way).
function along_files(root)
    out = Dict{String,String}()
    for f in readdir(joinpath(root, "csv"))
        m = match(r"^([A-Za-z]+)_[Pp]ara__([\d.,]+)%\.csv$", f)
        m === nothing && continue
        out["$(lowercase(m[1]))_$(replace(m[2], '.' => ','))%_para"] =
            joinpath(root, "csv", f)
    end
    return out
end

"""
    load_run(name; root) -> Run

Across the grain, parse one published file by its own header: the
thermocouple columns present (`2mm` … `40mm`, some missing), the mass in g
(`PM (g)`, written `PM (s)` in two files; absent in the dry runs), LF or CRLF
endings, `Temps` or `Time`. Along the grain (a name ending in `_para`), see
`load_along_run`.
"""
function load_run(name::AbstractString; root = nothing)
    endswith(name, "_para") &&
        return load_along_run(name; root = something(root, lemta_along_data()))
    root = something(root, lemta_data())
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
        :across,
        mc,
        rho,
        rep,
        t,
        depths,
        T,
        t,
        mass,
        cone_off(t, T[:, 1]),
        T_amb,
        T[1, :],
        defect,
    )
end

"""
    load_along_run(name; root) -> Run

One along-grain sheet: up to two blocks side by side, each opening with its
own `Time (h)` column, the thermocouples (`2mm` … `40mm`, some missing) at 1 s
and the mass (`M (g)`) on the balance's clock. The density is the deck's
along-grain one (Diarra et al. 2026 Table 1, //). Three mass records are
stamped 1 s apart where the balance logs every 1.07 s (`MASS_STEP`): on their
stamps poplar_21,3% loses its rate 10 min before its thermocouples show the
withdrawal, and the records span 7.0–7.1 h of a 7.5 h test. They are
restamped at `MASS_STEP` from the start of each stretch without a gap. The
withdrawal is read from the thermocouples, or from the mass where there are
none (`cone_off_mass`).
"""
function load_along_run(name::AbstractString; root = lemta_along_data())
    m = match(r"^([a-z]+)_([\d,]+)%_para$", name)
    m === nothing && error("unexpected along-grain run name $name")
    species = String(m[1])
    mc = parse(Float64, replace(m[2], ',' => '.')) / 100
    rho = Float64(deck(species)["along"]["rho_dry"])
    rows = [split(l, ',') for l in readlines(along_files(root)[name])]
    header = strip.(rows[1])
    starts = findall(==("Time (h)"), header)
    block(i0) =
        i0:(something(
            findnext(h->isempty(h)||h=="Time (h)", header, i0+1),
            length(header)+1,
        )-1)
    tcs = Tuple{Int,Int}[]
    itT = imass = itM = 0
    for i0 in starts
        cols = block(i0)
        k = findfirst(i -> header[i] == "M (g)", cols)
        if k === nothing
            itT = i0
            tcs = [
                (i, parse(Int, header[i][1:(end-2)])) for
                i in cols if occursin(r"^\d+mm$", header[i])
            ]
        else
            itM, imass = i0, cols[k]
        end
    end
    value(r, i) = i <= length(r) ? tryparse(Float64, r[i]) : nothing
    # one block's rows: those where every column of the block parses
    function read_block(it, cols)
        t, v = Float64[], Vector{Float64}[]
        for r in rows[2:end]
            x = [value(r, i) for i in (it, cols...)]
            any(isnothing, x) && continue
            push!(t, 3600 * x[1])
            push!(v, Float64.(x[2:end]))
        end
        return t, v
    end
    t_mass, mass = Float64[], Float64[]
    if imass > 0
        t_mass, v = read_block(itM, (imass,))
        mass = first.(v)
        t_mass = restamp(t_mass)
    end
    defect = get(DEFECTS, name, :none)
    if itT > 0
        t, v = read_block(itT, first.(tcs))
        T = reduce(vcat, permutedims.(v)) .+ 273.15
        depths = [d / 1000 for (_, d) in tcs]
        T_amb = mean(T[1:min(5, end), 1])
        T0 = T[1, :]
        t_off = cone_off(t, T[:, 1])
    else
        t, T, depths = t_mass, zeros(length(t_mass), 0), Float64[]
        T_amb, T0 = T_UNRECORDED, Float64[]
        t_off = defect === :mass_gap ? 9000.0 : cone_off_mass(t_mass, mass)
    end
    return Run(
        name,
        species,
        :along,
        mc,
        rho,
        0,
        t,
        depths,
        T,
        t_mass,
        mass,
        t_off,
        T_amb,
        T0,
        defect,
    )
end

# The balance's times restamped at MASS_STEP when the record steps 1 s, from
# the start of each stretch without a gap (a step over 5 s).
function restamp(t)
    length(t) < 2 && return t
    abs(median(diff(t)) - 1.0) < 0.01 || return t
    out = similar(t)
    s = 1
    for i in eachindex(t)
        i > 1 && t[i] - t[i-1] > 5 && (s = i)
        out[i] = t[s] + MASS_STEP * (t[i] - t[s])
    end
    return out
end

"""
The runs built from a deck alone, for the slice profiles that have no
recorded run: spruce along the grain at 30 % (Diarra et al. 2026 Fig. 11),
`T_UNRECORDED` throughout, under the cone for 2.5 h.
"""
slice_runs() = [
    Run(
        "spruce_30%_para_slices",
        "spruce",
        :along,
        0.30,
        Float64(deck("spruce")["along"]["rho_dry"]),
        0,
        collect(0.0:10.0:9000.0),
        Float64[],
        zeros(901, 0),
        Float64[],
        Float64[],
        Inf,
        T_UNRECORDED,
        Float64[],
        :none,
    ),
]

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

"""
The cone's withdrawal from the mass alone, for the runs without
thermocouples: the time in 2.2–3.2 h at which the loss rate, a straight line
fitted over the `window` after against the `window` before, rises most. On the
six along-grain runs that have both records it falls within 2 min of the
thermocouples' withdrawal.
"""
function cone_off_mass(t, m; window = 300.0)
    function slope(a, b)
        i, j = searchsortedfirst(t, a), searchsortedlast(t, b)
        j - i < 20 && return NaN
        x, y = view(t, i:j), view(m, i:j)
        x̄, ȳ = mean(x), mean(y)
        return sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
    end
    best, at = -Inf, Inf
    for x = 7920.0:30.0:11520.0
        jump = slope(x, x + window) - slope(x - window, x)
        isfinite(jump) && jump > best && ((best, at) = (jump, x))
    end
    return at
end

"""
The runs of one grain direction (`:across` or `:along`), the excluded one
dropped unless `excluded = true`.
"""
function load_runs(;
    grain::Symbol = :across,
    root = grain === :along ? lemta_along_data() : lemta_data(),
    excluded = false,
)
    runs = [load_run(n; root) for n in run_names(; grain, root)]
    excluded || filter!(r -> r.defect != :excluded, runs)
    return runs
end

"""
The model's initial water [g] in the 10 × 10 × 5 cm block. Across the grain,
MC/(1 + MC) × m₀ (dry basis, Diarra et al. Eq. 1) from the run's own first
mass, or MC ρ_dry V when the first mass is not the sample's. Along the grain
MC ρ_dry V, at the density `rho`: the first masses are 4–27 % below the
along-grain density over 500 cm³, so the block starts at the stated moisture
content and the mass is compared as the share of the water lost
([`water_fraction_lost`](@ref)).
"""
function initial_water_g(r::Run; rho = r.rho)
    r.mc == 0 && return 0.0
    (r.grain === :along || r.defect == :loss_only) && return r.mc * rho * VOLUME * 1000
    return r.mc / (1 + r.mc) * r.mass[1]
end

"""
The water the sample held when the test began [g]: MC/(1 + MC) × m₀, or
MC ρ_dry V when the first mass is not the sample's; NaN without a mass record.
"""
function measured_water_g(r::Run)
    r.mc == 0 && return 0.0
    r.defect == :loss_only && return r.mc * r.rho * VOLUME * 1000
    isempty(r.mass) && return NaN
    return r.mc / (1 + r.mc) * r.mass[1]
end

"""
Dry density the first mass gives over 500 cm³ [kg/m³], m₀/(1 + MC)/V; the
run's density without a mass record.
"""
mass_density(r::Run) = isempty(r.mass) ? r.rho : r.mass[1] / 1000 / (1 + r.mc) / VOLUME

# ---------------------------------------------------------------------------
# Decks and properties
# ---------------------------------------------------------------------------

deck(species::AbstractString) = TOML.parsefile(joinpath(DECKS, species * ".toml"))

"""
The deck of a run: the species' deck, along the grain with the entries of its
`[along]` table in place of the top level's (the density, the permeability,
the transfer coefficient `lambda` and the conductivity table `k_T_C`, `k`,
`rho_k`).
"""
deck(r::Run) =
    r.grain === :along ? merge(deck(r.species), deck(r.species)["along"]) : deck(r.species)

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
Conductivity along the grain [W/(m K)] of a deck holding the table `k_T_C`,
`k` (Flity et al. 2024 Table A.10, dry, measured on a sample of density
`rho_k`): the table at `T`, scaled by `ρ_dry/rho_k` (the conductivity is
linear in density at each temperature, their Fig. 9). Flity et al. publish no
law across species along the grain.
"""
k_along(d, ρ_dry, T) =
    interp(Float64.(d["k_T_C"]) .+ 273.15, Float64.(d["k"]), T) * ρ_dry /
    Float64(d["rho_k"])

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

"""
Saturation pressure of water [Pa] (Antoine, 1–374 °C), for the one-way laws'
diagnostic of where they evaporate into supersaturated vapour. The phase
change and the humid ambient use the model's own curve, `WATER_SATURATION`.
"""
function p_sat(T)
    θ = T - 273.15
    a, b, c = θ < 100 ? (8.07131, 1730.63, 233.426) : (8.14019, 1810.94, 244.485)
    return 133.322 * 10^(a - b / (c + θ))
end

"""
Heat taken by a one-way law (`LAWS`) desorbing one kg of bound water at `T`
[K] with `X` kg of water per kg of dry wood [J/kg]: the latent heat of free
water at `T` (linear, 2501 kJ/kg at 0 °C, 2257 at 100 °C) plus the
differential heat of sorption of examples/07, 875.6 exp(−11.22 X) kJ/kg,
scaled by `sorption`. Heats of reaction are quoted at the local temperature
(TRG §1.5.3), hence L(T). The phase change takes its heat from the
components' own heat capacities and its isotherm instead.
"""
heat_of_desorption(T, X; sorption = 1.0) =
    2.501e6 - 2369.0 * (T - 273.15) +
    sorption * 875.6441221735336e3 * exp(-11.21684273316275X)

"""
The one-way evaporation laws compared with the phase change, as (A [1/s],
E [J/mol], T_min [K], gate width [K]): `:arrhenius` is examples/07's set;
`:gated` the same switched on at 373.15 K over 2 K (kinetics audit §4.1a);
`:chan` the drying set the literature quotes from Chan, Kelbon & Krieger 1985
(kinetics audit §4.1). `:onset343` … `:onset363` and `:onset383` switch
examples/07's set on at 70, 80, 90 and 110 °C over 5 K, and `:onset373w10` at
100 °C over 10 K: a scan of where the latent sink sits, not candidate laws.
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

"""
The phase-change laws, the moisture model: `:sorption`, the bound water and
its vapour exchanging toward sorption equilibrium on the species' isotherm
(the deck's `isotherm`); `:freewater`, the same with the activity of free
water (`FreeLiquid`, no isotherm), which cannot hold water below saturation.
"""
const PHASE_CHANGE_LAWS = (:sorption, :freewater)

"""
The isotherms a deck can name (`isotherm`): the CAST sets of Dietenberger,
Glass & Boardman (2024) fitted to spruce and to beech, the first for the
softwoods and the second for the hardwoods.
"""
const ISOTHERMS = Dict("CAST_SPRUCE" => CAST_SPRUCE, "CAST_BEECH" => CAST_BEECH)

"""
The transfer coefficient and the ambient humidity at which the one-way laws
are compared: the dry wood's λ = 5e-7 m²/s and dry air. A one-way step makes
its vapour whatever the pore state and cannot take water back, so its
temperatures and mass hardly depend on either (permeability over two decades
moves the gated law by under 0.1 K and about 0.01 in mass).
"""
const ONE_WAY = (λ = 5e-7, ambient_rh = 0.0)

const IWOOD, IWATER, IVAPOUR, IAIR = 1, 2, 3, 4

"""
    material(d; law, κ, λ = 1.5e-6, sorption = 1.0, k_rule = PARALLEL,
             A_scale = 1.0, k_scale = 1.0, exchange_rate = 1e3,
             tortuosity = nothing)

The deck's wood. `law = :dry` gives the dry wood alone (thermal model only);
otherwise dry wood, bound water, water vapour and pore air. Under a law of
`PHASE_CHANGE_LAWS` the water and its vapour exchange two ways toward
equilibrium at the rate `exchange_rate` [1/s], with the wood the sorbing
solid; under a law of `LAWS` the water desorbs one way by that rate law.
The dry wood conducts by [`k_across`](@ref), or by [`k_along`](@ref) when the
deck holds the along-grain table (`deck(r)` of an along-grain run).
The dry wood's volume share is one at its dry density; the water holds its
own volume at 1000 kg/m³, and under `k_rule = PARALLEL` adds to the
conductivity by that share at 0.6 W/(m K), which gives the +8 to +26 % that
Flity et al. 2024 (Table 6) measured at about 10 % moisture across the grain;
`SERIES` removes the effect. `κ` [m²] and `λ` [m²/s] are the dry wood's Darcy
permeability and lumped vapour transfer coefficient; the water blocks neither
and carries no gas. A `tortuosity` replaces the lumped transfer by the
dusty-gas closure, the vapour diffusing through the pore air at
`D_VAPOUR_AIR` (scaled by T^1.75/P) reduced by φ/tortuosity. `A_scale`
multiplies a one-way law's pre-exponential factor and `sorption` the sorption
part of its heat; `k_scale` multiplies the dry wood's conductivity.
"""
function material(
    d;
    law::Symbol,
    κ::Real = 1e-15,
    λ::Real = 1.5e-6,
    sorption::Real = 1.0,
    k_rule = PARALLEL,
    A_scale::Real = 1.0,
    k_scale::Real = 1.0,
    exchange_rate::Real = 1e3,
    tortuosity::Union{Nothing,Real} = nothing,
)
    ρ = Float64(d["rho_dry"])
    Tε = Float64.(d["emissivity_T_C"]) .+ 273.15
    ε_tab = Float64.(d["emissivity"])
    εlaw = CallableProperty(T -> interp(Tε, ε_tab, T))
    k_dry = haskey(d, "k_T_C") ? (T -> k_along(d, ρ, T)) : (T -> k_across(ρ, T))
    wood = SolidComponent(
        :WOOD;
        ρ = ρ,
        c = CallableProperty(c_wood),
        k = CallableProperty(T -> k_scale * k_dry(T)),
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
    if law in PHASE_CHANGE_LAWS
        moisture = PhaseChange(
            :moisture,
            :MOISTURE => :WATER_VAPOR;
            activity = law === :freewater ? FreeLiquid() : ISOTHERMS[d["isotherm"]],
            sorbents = (:WOOD,),
            exchange_rate,
        )
        reactions, phase_changes = (), (moisture,)
    else
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
        reactions, phase_changes = (evaporation,), ()
    end
    return Material(
        name = Symbol(d["species"]),
        components = (wood, water, vapour, air),
        reactions = reactions,
        phase_changes = phase_changes,
        mixing = MaterialMixing(
            k = MixingSpec(k_rule),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
        gas_transport = tortuosity === nothing ? LumpedTransfer() :
                        DustyGas(;
            tortuosity,
            D_ref = D_VAPOUR_AIR,
            T_ref = 298.15,
            P_ref = P_ATM,
        ),
    )
end

"""
The pore gas at the start, (vapour, air) [kg/m³ bulk], over the condensed
concentrations `ξ0` at `T`: under a phase change the vapour in equilibrium
with the water and the air making up the ambient pressure, otherwise air
alone at the ambient pressure.
"""
function start_pore_gas(mat, ξ0, T)
    isempty(mat.phase_changes) &&
        return (0.0, background_gas_concentration(mat, :AIR, ξ0, T))
    p_eq = equilibrium_vapour_pressure(mat, ξ0, T)
    return (
        equilibrium_vapour_concentration(mat, ξ0, T),
        background_gas_concentration(mat, :AIR, ξ0, T; P = P_ATM - p_eq),
    )
end

# ---------------------------------------------------------------------------
# Scenario
# ---------------------------------------------------------------------------

# 1 under the cone, 0 after its withdrawal at `t_off`, a smooth step of width SWITCH.
cone_on(t, t_off) = isfinite(t_off) ? 0.5 * (1 - tanh((t - t_off) / SWITCH)) : 1.0

"""
    problem(r::Run; law = :sorption, cells = 100, t_end, flux = FLUX,
            rear_radiation = false, absorptivity, h_front, h_rear, κ, κ_scale,
            λ, ambient_rh = 0.5, exchange_rate = 1e3, tortuosity, rho,
            sorption, k_rule, A_scale, k_scale) -> (problem, mesh nodes)

The 1D slab of 50 mm, the grain across or along the heat flux as the run's:
the exposed face (`:top`, z = 50 mm) absorbs `absorptivity × flux` while the
cone is on and re-radiates at ε(T_s) to the ambient, with convection at the
deck's heating and cooling coefficients; the back (`:bottom`) loses heat at
the deck's `h_rear` (the whole loss Flity et al. 2025 fitted there) plus
radiation when `rear_radiation`. Both faces are open to the ambient air at the
relative humidity `ambient_rh`, its vapour mole fraction
`ambient_rh p_sat(T_a)/P` with `WATER_SATURATION`. The initial temperature is
the run's first row, interpolated in depth (uniform at `T_amb` without
thermocouples); under a phase change the pore vapour starts in equilibrium
with the water and the air makes up the ambient pressure. The deck is the
run's (`deck(r)`); `rho` replaces its dry density, the conductivity scaling
with it. `κ` defaults to the deck's `permeability`, scaled by `κ_scale`, and
`λ` to the deck's `lambda` (along the grain) or `LAMBDA`; the other keywords
are those of [`material`](@ref) or override the deck. With `calibrated`, the
incident flux and the conductivity take the factors fitted to the dry runs of
the species in the run's grain direction (calibration.toml, `<species>` and
`<species>_along`; without dry runs the deck stays).
"""
function problem(
    r::Run;
    law::Symbol = :sorption,
    cells::Int = 100,
    t_end::Real = r.t[end],
    flux::Real = FLUX,
    absorptivity::Union{Nothing,Real} = nothing,
    h_front::Union{Nothing,Tuple} = nothing,
    h_rear::Union{Nothing,Real} = nothing,
    rear_radiation::Bool = false,
    κ::Union{Nothing,Real} = nothing,
    κ_scale::Real = 1.0,
    λ::Union{Nothing,Real} = nothing,
    ambient_rh::Real = 0.5,
    exchange_rate::Real = 1e3,
    tortuosity::Union{Nothing,Real} = nothing,
    rho::Union{Nothing,Real} = nothing,
    sorption::Real = 1.0,
    k_rule = PARALLEL,
    A_scale::Real = 1.0,
    k_scale::Real = 1.0,
    calibrated::Bool = false,
    save_every::Real = 10.0,
)
    d = deck(r)
    rho === nothing || (d = merge(d, Dict("rho_dry" => Float64(rho))))
    ρ = Float64(d["rho_dry"])
    if calibrated
        c = get(
            TOML.parsefile(CALIBRATION),
            r.grain === :along ? r.species * "_along" : r.species,
            Dict("flux_factor" => 1.0, "k_factor" => 1.0),
        )
        flux *= c["flux_factor"]
        k_scale *= c["k_factor"]
    end
    κ = κ_scale * something(κ, Float64(get(d, "permeability", 1e-15)))
    λ = something(λ, Float64(get(d, "lambda", LAMBDA)))
    mat = material(
        d;
        law = r.mc == 0 ? :dry : law,
        κ,
        λ,
        sorption,
        k_rule,
        A_scale,
        k_scale,
        exchange_rate,
        tortuosity,
    )
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
    Tinit = isempty(dep) ? (z -> Ta) : (z -> interp(dep, T0, THICKNESS - z))

    if r.mc == 0
        bcs = BoundaryConditionSet(thermal = (top = top, bottom = bottom))
        ξinit = [z -> ρ]
    else
        x_v = ambient_rh * WATER_SATURATION(Ta) / P_ATM
        reservoir = GasReservoir(
            P = P_ATM,
            T = Ta,
            mole_fractions = (AIR = 1.0 - x_v, WATER_VAPOR = x_v),
        )
        bcs = BoundaryConditionSet(
            thermal = (top = top, bottom = bottom),
            mass = (top = AmbientGasBC(reservoir), bottom = AmbientGasBC(reservoir)),
        )
        ρw = initial_water_g(r; rho = ρ) / 1000 / VOLUME
        ξ0 = (ρ, ρw, 0.0, 0.0)
        ξinit = [
            z -> ξ0[1],
            z -> ξ0[2],
            z -> start_pore_gas(mat, ξ0, Tinit(z))[1],
            z -> start_pore_gas(mat, ξ0, Tinit(z))[2],
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
    simulate(r::Run; law = :sorption, kwargs...) -> NamedTuple

Solve the run's problem and return the series the comparison reads: the
thermocouple temperatures [K] (cell centres interpolated linearly to each
depth, the depths scaled by `probe_scale`), the sample mass [g] (dry wood and
water), the water left [g], the largest pore pressure and vapour-to-saturation
ratio, and the evaporation rate [g/s] with the share of it taking place where
the pore vapour pressure exceeds the saturation pressure at the cell
temperature (where a vapour-pressure-limited law would not evaporate).

Under a one-way law the evaporation is the reaction's rate. Under a phase
change it is the gross evaporation of the solution's extras, and
`phase_change` holds what the phase change adds: the gross condensation
[g/s], the largest ratio of the pore vapour pressure to its equilibrium
pressure where the wood holds water (X > 10⁻³), the largest moisture content,
the accepted and rejected steps, and the 10th, 50th and 90th percentiles of
the cell temperature weighted by the evaporation over the heating period [K].
`profiles` holds the moisture and pore state of every cell at 0, 0.5, 1, 1.5,
2 and 2.5 h and at the end (see [`write_profiles`](@ref)).
"""
function simulate(
    r::Run;
    law::Symbol = :sorption,
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
    two_way = moist && !isempty(mat.phase_changes)
    mass = zeros(nt)
    water = zeros(nt)
    pmax = fill(P_ATM, nt)
    smax = zeros(nt)
    evap = zeros(nt)
    evap_super = zeros(nt)
    rates = zeros(max(length(mat.reactions), 1))
    heats = similar(rates)
    cond = zeros(nt)
    ratio = zeros(nt)
    Xmax = zeros(nt)
    t_heat = min(9000.0, r.t[end], sol.t[end])
    T_evap, w_evap = Float64[], Float64[]
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
            if two_way
                e = max(cell_rate(mat, sol, ξ, c, k), 0.0) * widths[c] * AREA * 1000
                X = Pyrolysis.Physics.moisture_content(ξ[IWATER], max(ξ[IWOOD], 0.0))
                Xmax[k] = max(Xmax[k], X)
                p_v = sol.gas.partial_pressure[vapour_gas(sol), c, k]
                X > 1e-3 &&
                    (ratio[k] = max(ratio[k], p_v / equilibrium_vapour_pressure(mat, ξ, T)))
                if e > 0 && sol.t[k] <= t_heat
                    push!(T_evap, T)
                    push!(w_evap, e)
                end
            else
                Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, mat)
                e = rates[1] * widths[c] * AREA * 1000
                evap[k] += e
            end
            s > 1 && (evap_super[k] += e)
        end
    end
    if two_way
        evap .= sol.extras.evaporation_rate .* AREA .* 1000
        cond .= sol.extras.condensation_rate .* AREA .* 1000
    end
    phase_change =
        two_way ?
        (;
            condensation = cond,
            pv_over_peq = ratio,
            X_max = Xmax,
            accepted = sol.statistics.accepted,
            rejected = sol.statistics.rejected,
            T_evap = map(q -> weighted_quantile(T_evap, w_evap, q), (0.1, 0.5, 0.9)),
        ) : nothing
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
        phase_change,
        profiles = moist ? profiles(r, sol, nodes) : nothing,
    )
end

# Index of the water vapour among the solution's gas species.
vapour_gas(sol) = findfirst(==(:WATER_VAPOR), sol.gas.species)

# The phase change's rate [kg/(m³ s)] in cell `c` at saved time `k`, read as
# the kernel reads it: the signed concentrations `ξ`, the sorbent at its
# positive part, the cell's pore fraction.
function cell_rate(mat, sol, ξ, c, k)
    pc = only(mat.phase_changes)
    X = Pyrolysis.Physics.moisture_content(pc, ξ, map(x -> max(x, 0.0), ξ))
    return Pyrolysis.Physics.phase_change_rate(
        pc,
        X,
        sol.T[c, k],
        sol.gas.porosity[c, k],
        ξ[pc.vapour],
    )
end

"""The `q` quantile of the values `x` weighted by `w` (NaN when empty)."""
function weighted_quantile(x, w, q)
    isempty(x) && return NaN
    p = sortperm(x)
    cs = cumsum(w[p])
    return x[p][findfirst(>=(q * cs[end]), cs)]
end

# The profiles of the moisture content and the pore state, cell by cell at
# the times `PROFILE_HOURS` and at the end.
const PROFILE_HOURS = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5)

function profiles(r::Run, sol, nodes)
    mat = sol.problem.material
    NC = size(sol.ξ, 1)
    nt = length(sol.t)
    law =
        isempty(mat.phase_changes) ? ISOTHERMS[deck(r.species)["isotherm"]] :
        only(mat.phase_changes).activity
    rows = NTuple{7,Float64}[]
    for th in (PROFILE_HOURS..., sol.t[end] / 3600)
        k = something(findfirst(>=(3600th - 1e-6), sol.t), nt)
        for c in reverse(eachindex(diff(nodes)))
            ξ = ntuple(j -> sol.ξ[j, c, k], NC)
            T = sol.T[c, k]
            X = Pyrolysis.Physics.moisture_content(ξ[IWATER], max(ξ[IWOOD], 0.0))
            push!(
                rows,
                (
                    sol.t[k] / 3600,
                    1000 * (nodes[end] - 0.5 * (nodes[c] + nodes[c+1])),
                    T,
                    X,
                    sol.gas.partial_pressure[vapour_gas(sol), c, k],
                    Pyrolysis.Physics.equilibrium_vapour_pressure(
                        law,
                        WATER_SATURATION,
                        X,
                        T,
                    ),
                    sol.gas.pressure[c, k],
                ),
            )
        end
    end
    return rows
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
# Moist runs under the phase change: the deck as published (`sorption`) and
# with the dry runs' factors (`sorption_cal`), and on the calibrated deck the
# transfer coefficient's ladder, the ambient humidity, the vent permeability,
# the exchange rate, the numerics, the dusty-gas closure and free water.
VARIANTS["sorption"] = (runs = :moist, kwargs = r -> (; law = :sorption))
VARIANTS["sorption_cal"] =
    (runs = :moist, kwargs = r -> (; law = :sorption, calibrated = true))
for (tag, kw) in (
    "lambda5e-7" => (λ = 5e-7,),
    "lambda5e-6" => (λ = 5e-6,),
    "rh0" => (ambient_rh = 0.0,),
    "rh65" => (ambient_rh = 0.65,),
    "vent" => (κ = 1e-15,),
    "K1e2" => (exchange_rate = 1e2,),
    "K1e4" => (exchange_rate = 1e4,),
    "cells200" => (cells = 200,),
    "tight" => (reltol = 1e-8, dtmax = 10.0),
    "dusty4" => (tortuosity = 4.0,),
    "dusty12" => (tortuosity = 12.0,),
    "freewater" => (law = :freewater,),
)
    VARIANTS["sorption_cal_"*tag] =
        (runs = :moist, kwargs = r -> (; law = :sorption, calibrated = true, kw...))
end
# Along the grain: the transfer coefficient at the across-grain value and at
# 5e-5 (with the ladder's 5e-6, a third to three times the along-grain 1.5e-5),
# the permeability across the grain in place of the along-grain one, and the
# dry density the first mass gives (`mass_density`) in place of Table 1's.
for (tag, kw) in (
    "lambda1.5e-6" => r -> (λ = 1.5e-6,),
    "lambda5e-5" => r -> (λ = 5e-5,),
    "kappa_across" => r -> (κ = Float64(deck(r.species)["permeability"]),),
    "rho_mass" => r -> (rho = mass_density(r),),
)
    VARIANTS["sorption_cal_"*tag] =
        (runs = :moist, kwargs = r -> (; law = :sorption, calibrated = true, kw(r)...))
end
# The one-way laws (at the settings of `ONE_WAY`): the deck as published
# (`<law>`) and with the dry runs' factors (`<law>_cal`); the sensitivities on
# the calibrated deck.
for law in (:arrhenius, :gated, :chan)
    l = String(law)
    VARIANTS[l] = (runs = :moist, kwargs = r -> (; law, ONE_WAY...))
    VARIANTS[l*"_cal"] =
        (runs = :moist, kwargs = r -> (; law, calibrated = true, ONE_WAY...))
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
            (runs = :moist, kwargs = r -> (; law, calibrated = true, ONE_WAY..., kw...))
    end
end
for law in (:onset343, :onset353, :onset363, :onset373w10, :onset383)
    VARIANTS[String(law)*"_cal"] =
        (runs = :moist, kwargs = r -> (; law, calibrated = true, ONE_WAY...))
end

# ---------------------------------------------------------------------------
# Series files
# ---------------------------------------------------------------------------

series_path(variant, name) = joinpath(RESULTS, variant, name * ".csv")
profiles_path(variant, name) = joinpath(RESULTS, variant, name * "_profiles.csv")

"""
    write_series(path, s)

The series of [`simulate`](@ref), one row per saved time. A phase change adds
the columns `condensation_g_s`, `pv_over_peq_max` and `X_max` and a line
`# phase_change accepted=… rejected=… T_evap_p10_K=… T_evap_p50_K=…
T_evap_p90_K=…` after the line `# retcode=… solve_time=…`.
"""
function write_series(path, s)
    mkpath(dirname(path))
    pc = s.phase_change
    open(path, "w") do io
        tcs = [@sprintf("T_%gmm_K", 1000d) for d in s.depths]
        extra = pc === nothing ? String[] : ["condensation_g_s", "pv_over_peq_max", "X_max"]
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
                    ];
                    extra
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
                ];
                pc === nothing ? Float64[] :
                [pc.condensation[k], pc.pv_over_peq[k], pc.X_max[k]]
            ]
            println(io, join((@sprintf("%.6g", v) for v in vals), ","))
        end
        println(io, "# retcode=$(s.retcode) solve_time=$(round(s.solve_time, digits = 2))")
        pc === nothing || @printf(
            io,
            "# phase_change accepted=%d rejected=%d T_evap_p10_K=%.2f T_evap_p50_K=%.2f T_evap_p90_K=%.2f\n",
            pc.accepted,
            pc.rejected,
            pc.T_evap...
        )
    end
    return path
end

"""
    read_series(path) -> NamedTuple

A series written by [`write_series`](@ref). `phase_change` is `nothing` for a
one-way law, and otherwise holds the phase change's columns (`condensation`,
`pv_over_peq`, `X_max`) and its line's values (`accepted`, `rejected`,
`T_evap`, the three percentiles).
"""
function read_series(path)
    all_lines = readlines(path)
    lines = filter(l -> !startswith(l, "#"), all_lines)
    header = split(lines[1], ',')
    data = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
    col(name) = [row[findfirst(==(name), header)] for row in data]
    tcs = [h for h in header if occursin(r"^T_[\d.]+mm_K$", h)]
    depths = [parse(Float64, h[3:(end-4)]) / 1000 for h in tcs]
    T = isempty(tcs) ? zeros(length(data), 0) : reduce(hcat, [col(h) for h in tcs])
    meta = filter(l -> startswith(l, "# retcode"), all_lines)
    pcline = filter(l -> startswith(l, "# phase_change"), all_lines)
    phase_change = if isempty(pcline)
        nothing
    else
        v = Dict(
            String(a) => parse(Float64, b) for
            (a, b) in (split(f, '=') for f in split(pcline[1])[3:end])
        )
        (;
            condensation = col("condensation_g_s"),
            pv_over_peq = col("pv_over_peq_max"),
            X_max = col("X_max"),
            accepted = round(Int, v["accepted"]),
            rejected = round(Int, v["rejected"]),
            T_evap = (v["T_evap_p10_K"], v["T_evap_p50_K"], v["T_evap_p90_K"]),
        )
    end
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
        phase_change,
    )
end

"""
    write_profiles(path, s)

The profiles of [`simulate`](@ref), one row per cell and time: the time [h],
the cell centre's depth below the heated face [mm], T [K], the moisture
content X (water over dry wood), the vapour's partial pressure, its
equilibrium pressure on the run's activity law (the species' isotherm under a
one-way law) and the pore pressure [Pa].
"""
function write_profiles(path, s)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "t_h,depth_mm,T_K,X,p_v_Pa,p_eq_Pa,P_Pa")
        for row in s.profiles
            @printf(io, "%.3f,%.2f,%.2f,%.5f,%.1f,%.1f,%.1f\n", row...)
        end
    end
    return path
end

"""Profiles written by [`write_profiles`](@ref): Dict(time [h] => rows (depth_mm, X))."""
function read_profiles(path)
    out = Dict{Float64,Vector{Tuple{Float64,Float64}}}()
    for l in readlines(path)[2:end]
        f = parse.(Float64, split(l, ','))
        push!(get!(out, round(f[1], digits = 2), Tuple{Float64,Float64}[]), (f[2], f[4]))
    end
    return out
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
of the record within 30 s on the balance's clock (the balance has
single-sample spikes of up to 2 g; NaN where it recorded nothing), or in a
model series.
"""
function mass_loss(r::Run, t)
    isempty(r.mass) && return NaN
    function window(t0)
        k = findall(x -> abs(x - t0) <= 30, r.t_mass)
        return isempty(k) ? NaN : median(r.mass[k])
    end
    return window(0.0) - window(t)
end
mass_loss(s::NamedTuple, t) = s.mass[1] - interp(s.t, s.mass, t)

"""
Share of the initial water lost by time `t`: the mass lost over
[`measured_water_g`](@ref) in the measurement, over the water the model
started with in a model series.
"""
water_fraction_lost(r::Run, t) = mass_loss(r, t) / measured_water_g(r)
water_fraction_lost(s::NamedTuple, t) = mass_loss(s, t) / s.water[1]

"""
Thermocouple accuracy [K] (the IEC 60584 class Diarra et al. quote): ±2 K up
to 100 °C, ±1.5 % of the reading in °C above.
"""
tc_accuracy(T) = max(2.0, 0.015 * (T - 273.15))

"""
Temperature errors of a run against a series over `window`: per depth [mm]
(rms K, bias K, NRMS), all depths pooled (`rms`), and the share of the
readings the model meets within the thermocouple accuracy (`within`).
"""
function temperature_errors(r::Run, s; window = (0.0, Inf))
    t_hi = min(window[2], s.t[end], r.t[end])
    keep = findall(t -> window[1] <= t <= t_hi, r.t)
    per = Dict{Int,NamedTuple}()
    pooled = Float64[]
    within = 0
    for (j, d) in enumerate(r.depths)
        js = findfirst(x -> isapprox(x, d; atol = 1e-6), s.depths)
        y = r.T[keep, j]
        m = [interp(s.t, view(s.T, :, js), t) for t in r.t[keep]]
        e = m .- y
        append!(pooled, e)
        within += count(abs.(e) .<= tc_accuracy.(y))
        per[round(Int, 1000d)] = (
            rms = sqrt(mean(e .^ 2)),
            bias = mean(e),
            nrms = sqrt(mean(e .^ 2)) / (maximum(y) - minimum(y)),
        )
    end
    isempty(pooled) && return (; per, rms = NaN, within = NaN)
    return (; per, rms = sqrt(mean(pooled .^ 2)), within = within / length(pooled))
end

"""Trapezoidal integral of `y` over `t`."""
trapz(t, y) = sum(0.5 * (y[i] + y[i+1]) * (t[i+1] - t[i]) for i = 1:(length(t)-1))

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
