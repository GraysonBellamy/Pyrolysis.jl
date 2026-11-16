# The LEMTA dry argon case (doi:10.57745/QRYAJ9): oven-dry wood blocks of seven
# species in an argon-filled chamber under a cone heater at 25, 35 and 45 kW/m²,
# with the mass of one block per species and flux, and in-depth temperatures of
# another at 35 and 45 kW/m². Data loading and the defect list, the species
# decks to Pyrolysis.jl problems, and the model outputs at the thermocouple
# depths. See README.md.
module LemtaArgon

using Pyrolysis
using Printf
using SHA
using TOML

export lemta_argon_data,
    check_data,
    MassRun,
    TemperatureRun,
    load_mass,
    load_temperature,
    mass_runs,
    temperature_runs,
    usable,
    live,
    SPECIES,
    DEPTHS,
    RESULTS,
    chamber,
    deck,
    material,
    problem,
    simulate,
    char_yield,
    reaction_heats,
    interp,
    silcal_data

const HERE = @__DIR__
const DECKS = joinpath(HERE, "decks")
const RESULTS = joinpath(HERE, "results")

const VOLUME = 0.1 * 0.1 * 0.05   # m³, the 10 × 10 × 5 cm block
const DEPTHS = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40] ./ 1000   # m, nominal

# The dataset's species names and the decks'.
const SPECIES = Dict(
    "Balsa" => "balsa",
    "Poplar" => "poplar",
    "Fir" => "fir",
    "Larch" => "larch",
    "Oak" => "oak",
    "Beech" => "beech",
    "Ipe" => "ipe",
)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

lemta_argon_data() =
    get(ENV, "LEMTA_ARGON_DATA", joinpath(homedir(), "Documents", "data", "lemta-QRYAJ9"))

"""
Check every file of the dataset against its pinned SHA-256 (data_files.toml):
each original upload, and each spreadsheet's tab-delimited export.
"""
function check_data(; root = lemta_argon_data())
    isdir(root) || error(
        "LEMTA argon data not found at $root (set LEMTA_ARGON_DATA, or run fetch_data.sh)",
    )
    for f in TOML.parsefile(joinpath(HERE, "data_files.toml"))["files"]
        for (name, sum) in
            ((f["name"], f["sha256"]), (get(f, "tab", ""), get(f, "tab_sha256", "")))
            isempty(name) && continue
            path = joinpath(root, name)
            isfile(path) || error("missing $name in $root (run fetch_data.sh)")
            got = bytes2hex(open(sha256, path))
            got == sum || error("$name differs from the published file (SHA-256 $got)")
        end
    end
    return root
end

"""
Files whose name gives the wrong heat flux: each matches, curve for curve, the
thesis's figure of the other flux (Flity 2023, Annex E).
- `Mass_Balsa_35kW.m-2_143kg.m-3` is Fig. E.1(e, f), 45 kW/m²: m/m₀ 0.645 and
  3.3 g/(m² s) at 600 s, and the same dip of the rate at the start;
  `Mass_Balsa_35kW.m-2_76kg.m-3` is Fig. E.1(c, d). The set has no balsa mass
  file named 45.
- `Temperature_Ipe_35kW.m-2_842kg.m-3` is Fig. E.8(g), 45 kW/m²: the sample
  starts at 48 °C throughout, and its 6 and 10 mm records cross at 497–505 s;
  `Temperature_Ipe_35kW.m-2_885kg.m-3` is Fig. 5.13(g). The set has no ipé
  temperature file named 45.
"""
const FLUX = Dict(
    "Mass_Balsa_35kW.m-2_143kg.m-3" => 45.0,
    "Temperature_Ipe_35kW.m-2_842kg.m-3" => 45.0,
)

# Mass_<Species>_<flux>kW.m-2_<density>kg.m-3 and Temperature_…: the species,
# the heat flux [kW/m²] (corrected by `FLUX`) and the oven-dry density.
function parse_name(name)
    m = match(r"^(Mass|Temperature)_([A-Za-z]+)_([\d.]+)kW\.m-2_([\d.]+)kg\.m-3$", name)
    m === nothing && error("unexpected file name $name")
    haskey(SPECIES, m.captures[2]) || error("$name: unknown species $(m.captures[2])")
    flux = get(FLUX, name, parse(Float64, m.captures[3]))
    return SPECIES[m.captures[2]], flux, parse(Float64, m.captures[4])
end

# A header cell with spaces removed and lower-cased: the export drops the
# spreadsheet's spaces ("Temps (s) " is "Temps(s)").
_key(h) = lowercase(replace(h, r"\s" => ""))

function read_tab(path)
    lines = readlines(path)
    header = split(chomp(replace(lines[1], '\r' => "")), '\t')
    rows = Vector{Vector{Float64}}()
    for l in lines[2:end]
        f = split(chomp(replace(l, '\r' => "")), '\t')
        length(f) == length(header) ||
            error("$path: a row of $(length(f)) fields under $(length(header)) headers")
        push!(rows, [isempty(strip(x)) ? NaN : parse(Float64, x) for x in f])
    end
    return header, reduce(vcat, permutedims.(rows))
end

"""
One mass test as published: `m` is the mass over the initial mass and `mlr`
the mass loss rate per unit area [g/(m² s)], both as published, against `t`
[s] from the start of exposure. `rho` is the oven-dry density of the file
name and `m0` = ρ V [g], which the published rate integrates back to within
1.5 %. `flux` is the heat flux [kW/m²] the test was run at (see `FLUX`).
`t_settle` is the end of the balance's disturbance at the start: the last
time m/m₀ exceeds one (zero if it never does), up to 29 s in these files.
"""
struct MassRun
    name::String
    species::String
    flux::Float64
    rho::Float64
    m0::Float64
    t::Vector{Float64}
    m::Vector{Float64}
    mlr::Vector{Float64}
    t_settle::Float64
end

Base.show(io::IO, r::MassRun) = print(
    io,
    "MassRun($(r.name): $(r.species), $(r.flux) kW/m², ρ $(r.rho), ",
    "$(round(r.t[end], digits = 0)) s, m/m₀ end $(round(r.m[end], digits = 3)), settled $(r.t_settle) s)",
)

"""
    load_mass(name; root) -> MassRun

Parse one mass file by its own header: a time column (`Temps (s)`), the
normalised mass (`Perte en masse normalisee (-)`, which is m/m₀, not the loss)
and the rate (`MLRPUA (g/(m2.s))`).
"""
function load_mass(name::AbstractString; root = lemta_argon_data())
    species, flux, rho = parse_name(name)
    header, A = read_tab(joinpath(root, name * ".tab"))
    k = _key.(header)
    it = findfirst(h -> startswith(h, "temps(s)") || startswith(h, "time(s)"), k)
    im = findfirst(h -> occursin("normalis", h), k)
    iq = findfirst(h -> startswith(h, "mlrpua"), k)
    any(isnothing, (it, im, iq)) && error("$name: unexpected header $header")
    any(isnan, A) && error("$name: empty cells")
    t, m = A[:, it], A[:, im]
    k1 = findlast(>(1.0), m)
    t_settle = k1 === nothing ? 0.0 : t[k1]
    return MassRun(name, species, flux, rho, rho * VOLUME * 1000, t, m, A[:, iq], t_settle)
end

"""
One in-depth temperature test as published. `T` [K] has one column per
nominal depth in `DEPTHS` (`NaN` where a thermocouple gave no reading), against
`t` [s] from the start of exposure. `T0` is the first row, the initial profile:
several samples start warm and graded (balsa at 45 kW/m²: 68 °C at 2 mm, 37 °C
at 40 mm). `t_end` is the end of exposure: the record's end unless the record
shows the heater withdrawn. `cuts[j]` lists the time windows in which the
thermocouple at `DEPTHS[j]` is not used (`(-Inf, Inf)` for a dead one), and
`notes` says why: see `DEFECTS` and `GLUE_FAILED`.
"""
struct TemperatureRun
    name::String
    species::String
    flux::Float64
    rho::Float64
    t::Vector{Float64}
    T::Matrix{Float64}
    T0::Vector{Float64}
    t_end::Float64
    cuts::Vector{Vector{Tuple{Float64,Float64}}}
    notes::Vector{String}
end

live(r::TemperatureRun, j::Int) = !any(c -> c == (-Inf, Inf), r.cuts[j])

Base.show(io::IO, r::TemperatureRun) = print(
    io,
    "TemperatureRun($(r.name): $(r.species), $(r.flux) kW/m², ρ $(r.rho), ",
    "$(count(j -> live(r, j), eachindex(DEPTHS))) of 12 TCs, $(round(r.t[end], digits = 0)) s)",
)

"""
Defects of the published temperature records found by reading every file
(check_data.jl draws them): `(depth_mm, reason, (t_from, t_to))`, a window [s]
in which that thermocouple is not used, or `(:exposure_end, reason, (t, Inf))`
when the record goes on after the heater was withdrawn. Columns without a
single reading are found by the loader and are not listed. The crosstalk
windows are steps of equal size and opposite sign on two channels at the same
instants.
"""
const DEFECTS = Dict{String,Vector{Tuple{Any,String,Tuple{Float64,Float64}}}}(
    "Temperature_Balsa_45kW.m-2_102kg.m-3" => [
        (6, "crosstalk with 10 mm", (400.0, 412.0)),
        (10, "crosstalk with 6 mm", (400.0, 412.0)),
    ],
    "Temperature_Beech_35kW.m-2_600kg.m-3" => [
        (2, "crosstalk with 6 mm", (965.0, 1005.0)),
        (6, "crosstalk with 2 mm", (965.0, 1005.0)),
        (
            :exposure_end,
            "every thermocouple to 10 mm falls 12–22 K in 2 s at 1079 s",
            (1078.0, Inf),
        ),
    ],
    "Temperature_Ipe_35kW.m-2_842kg.m-3" => [
        (6, "crosstalk with 10 mm", (470.0, 510.0)),
        (10, "crosstalk with 6 mm", (470.0, 510.0)),
    ],
    "Temperature_Larch_45kW.m-2_630kg.m-3" => [
        (2, "crosstalk with 6 mm from 562 s to the end", (560.0, Inf)),
        (6, "crosstalk with 2 mm from 562 s to the end", (560.0, Inf)),
    ],
)

"""
Species whose glue line opened in the cone tests, so that the heater shone on
the thermocouples near the face, which then read above the wood's temperature
(Flity 2023 thesis §5.4, printed p.110: balsa, spruce, oak and above all
beech; poplar, larch and ipé held). Their 2, 4 and 6 mm records are not used
above `GLUE_LIMIT`.
"""
const GLUE_FAILED = ("balsa", "oak", "beech")
const GLUE_LIMIT = 673.15   # K

"""
    load_temperature(name; root) -> TemperatureRun

Parse one temperature file by its own header: a time column (`Temps (s)`) and
one column per thermocouple, named by its depth (`2mm` … `40mm`). Blank cells
are `NaN`. A column with no reading is dead; `DEFECTS` and `GLUE_FAILED` cut
the others.
"""
function load_temperature(name::AbstractString; root = lemta_argon_data())
    species, flux, rho = parse_name(name)
    header, A = read_tab(joinpath(root, name * ".tab"))
    k = _key.(header)
    it = findfirst(h -> startswith(h, "temps(s)") || startswith(h, "time(s)"), k)
    it === nothing && error("$name: no time column in $header")
    t = A[:, it]
    T = fill(NaN, length(t), length(DEPTHS))
    for (i, h) in enumerate(k)
        m = match(r"^(\d+)mm$", h)
        m === nothing && continue
        j = findfirst(d -> isapprox(d, parse(Int, m.captures[1]) / 1000), DEPTHS)
        j === nothing && error("$name: thermocouple at $h is not a nominal depth")
        T[:, j] = A[:, i] .+ 273.15
    end
    cuts = [Tuple{Float64,Float64}[] for _ in DEPTHS]
    notes = String[]
    for j in eachindex(DEPTHS)
        mm = round(Int, 1000DEPTHS[j])
        if all(isnan, view(T, :, j))
            push!(cuts[j], (-Inf, Inf))
            push!(notes, "$mm mm: no reading")
        elseif species in GLUE_FAILED && mm <= 6
            i = findfirst(>=(GLUE_LIMIT), view(T, :, j))
            if i !== nothing
                push!(cuts[j], (t[i], Inf))
                push!(
                    notes,
                    @sprintf("%d mm: glue line, cut above 400 °C (from %g s)", mm, t[i])
                )
            end
        end
    end
    t_end = t[end]
    for (what, reason, window) in get(DEFECTS, name, [])
        if what === :exposure_end
            t_end = window[1]
        else
            push!(cuts[findfirst(d -> isapprox(d, what / 1000), DEPTHS)], window)
        end
        push!(
            notes,
            what === :exposure_end ? "exposure ends: $reason" : "$what mm: $reason",
        )
    end
    return TemperatureRun(name, species, flux, rho, t, T, T[1, :], t_end, cuts, notes)
end

"""The temperature at depth `j`, `NaN` where the thermocouple is not usable."""
function usable(r::TemperatureRun, j::Int)
    v = r.T[:, j]
    for (a, b) in r.cuts[j]
        v[a .<= r.t .<= b] .= NaN
    end
    v[r.t .> r.t_end] .= NaN
    return v
end

_names(prefix; root) = sort([
    splitext(f)[1] for f in readdir(root) if startswith(f, prefix) && endswith(f, ".tab")
])

mass_runs(; root = lemta_argon_data()) = [load_mass(n; root) for n in _names("Mass_"; root)]
temperature_runs(; root = lemta_argon_data()) =
    [load_temperature(n; root) for n in _names("Temperature_"; root)]

# ---------------------------------------------------------------------------
# Decks and properties
# ---------------------------------------------------------------------------

const P_ATM = 101325.0          # Pa
const σ_SB = 5.670374419e-8     # W/(m² K⁴)

"""The chamber and the laws every species shares (decks/chamber.toml)."""
chamber() = TOML.parsefile(joinpath(DECKS, "chamber.toml"))

"""A species deck (decks/<species>.toml)."""
deck(species::AbstractString) = TOML.parsefile(joinpath(DECKS, species * ".toml"))

"""Piecewise linear in `xs`, held constant outside them."""
function interp(xs::AbstractVector, ys::AbstractVector, x::Real)
    x <= xs[1] && return float(ys[1])
    x >= xs[end] && return float(ys[end])
    i = searchsortedlast(xs, x)
    w = (x - xs[i]) / (xs[i+1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i+1]
end

"""Final char yield of the kinetics, α2γ2 + α3γ3 (kg per kg of dry wood)."""
function char_yield(d)
    α, γ = d["kinetics"]["alpha"], d["kinetics"]["gamma"]
    return α[2] * γ[2] + α[3] * γ[3]
end

"""
Heats of the four reactions per kg of their reactant [J/kg, positive
endothermic], from the deck's heats per g of initial dry wood: H1/α1, H2/α2,
0 and H4/(α2γ2).
"""
function reaction_heats(d)
    α, γ = d["kinetics"]["alpha"], d["kinetics"]["gamma"]
    h = d["heats"]
    return (1e3h["H1"] / α[1], 1e3h["H2"] / α[2], 0.0, 1e3h["H4"] / (α[2] * γ[2]))
end

# Linear law of the temperature in °C, with T in K.
_celsius_law(c) = T -> c[1] + c[2] * (T - 273.15)

"""Dry wood across the grain [W/(m K)] at dry density `ρ` and `T` [K]."""
k_wood(ch, ρ, T) = ch["k_wood"][1] + ch["k_wood"][2] * (ρ / 1000) * sqrt(T / 433)

"""Char across the grain [W/(m K)] at `T` [K], by the named law of [k_char]."""
function k_char(ch, law, T)
    p = ch["k_char"][law]
    return p["a"] * T + p["b"] + 16σ_SB * T^3 / (3p["beta"])
end

"""
    material(d, rho; ch = chamber(), k_rule = PARALLEL, volatiles_c, permeability)

The species `d` at dry density `rho` [kg/m³]: the three pseudo-components of
the dry wood, the chars C2, C3 and the heat-only C4, the volatiles and the
chamber's argon.
- Each pseudo-component has the wood's bulk density, so the three fill the
  volume in their shares α. Each char has its reactant's volume (C2 and C4 at
  γ2 ρ, C3 at γ3 ρ), on a fixed mesh: the solid does not shrink, as in the
  thesis model.
- The wood components carry the wood's c, k(ρ, T), ε(T) and absorptivity; the
  chars the char's c and k law and ε = α = 0.95. Conductivity mixes by
  `k_rule` over the condensed volume shares, emissivity and absorptivity by
  those shares.
- Every solid has the chamber deck's permeability and lumped transfer
  coefficient; the volatiles' heat capacity is `volatiles_c`.
"""
function material(
    d,
    rho::Real;
    ch = chamber(),
    k_rule = PARALLEL,
    volatiles_c::Real = ch["volatiles_c"],
    permeability::Real = ch["permeability"],
)
    kin = d["kinetics"]
    α, γ = kin["alpha"], kin["gamma"]
    Tε = Float64.(d["emissivity_T_C"]) .+ 273.15
    εtab = Float64.(d["emissivity"])
    c_w = CallableProperty(_celsius_law(ch["c_wood"]))
    c_c = CallableProperty(_celsius_law(ch["c_char"]))
    k_w = CallableProperty(T -> k_wood(ch, rho, T))
    law = d["char_conductivity"]
    k_c = CallableProperty(T -> k_char(ch, law, T))
    ε_w = CallableProperty(T -> interp(Tε, εtab, T))
    ρs = ch["skeletal_density"]
    wood(name) = SolidComponent(
        name;
        ρ = rho,
        c = c_w,
        k = k_w,
        ε = ε_w,
        absorptivity = d["absorptivity"],
        λ = ch["transfer_wood"],
        κ = permeability,
        ρ_intrinsic = ρs,
    )
    char(name, yield) = SolidComponent(
        name;
        ρ = yield * rho,
        c = c_c,
        k = k_c,
        ε = ch["char_emissivity"],
        absorptivity = ch["char_absorptivity"],
        λ = ch["transfer_char"],
        κ = permeability,
        ρ_intrinsic = ρs,
    )
    components = (
        wood(:WOOD_1),
        wood(:WOOD_2),
        wood(:WOOD_3),
        char(:CHAR_2, γ[2]),
        char(:CHAR_3, γ[3]),
        char(:CHAR_4, γ[2]),
        GaseousComponent(
            :VOLATILES;
            M = ch["volatiles_M"],
            c = volatiles_c,
            k = ch["volatiles_k"],
        ),
        GaseousComponent(:ARGON; M = ch["argon_M"], c = ch["argon_c"], k = ch["argon_k"]),
    )
    A, E = kin["A"], kin["E"]
    H = reaction_heats(d)
    reactions = (
        Reaction(:R1, :WOOD_1 => :VOLATILES; A = A[1], E = E[1], h = H[1]),
        Reaction(
            :R2,
            :WOOD_2 => (:CHAR_2, :VOLATILES);
            A = A[2],
            E = E[2],
            h = H[2],
            yields = (γ[2], 1 - γ[2]),
        ),
        Reaction(
            :R3,
            :WOOD_3 => (:CHAR_3, :VOLATILES);
            A = A[3],
            E = E[3],
            h = H[3],
            yields = (γ[3], 1 - γ[3]),
        ),
        Reaction(:R4, :CHAR_2 => :CHAR_4; A = ch["A4"], E = ch["E4"], h = H[4]),
    )
    return Material(
        name = Symbol(d["species"]),
        components = components,
        reactions = reactions,
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

"""
The calcium-silicate runs in the chamber (decks/silcal.csv, Th Fig. 5.5): by
nominal flux [kW/m²], `(t, depths, T)` with `t` [s], `depths` [m] from the
exposed face (0.05 is the rear face) and `T` [K], `NaN` before a curve starts.
"""
function silcal_data()
    lines = filter(l -> !startswith(l, "#"), readlines(joinpath(DECKS, "silcal.csv")))
    header = split(lines[1], ',')
    depths =
        [parse(Int, match(r"T_(\d+)mm_C", h).captures[1]) / 1000 for h in header[3:end]]
    out = Dict{Int,Any}()
    for flux in unique(parse(Int, split(l, ',')[1]) for l in lines[2:end])
        rows = [split(l, ',') for l in lines[2:end] if parse(Int, split(l, ',')[1]) == flux]
        t = [parse(Float64, r[2]) for r in rows]
        T = reduce(
            vcat,
            [
                permutedims([
                    isempty(x) ? NaN : parse(Float64, x) + 273.15 for x in r[3:end]
                ]) for r in rows
            ],
        )
        out[flux] = (; t, depths, T)
    end
    return out
end

"""
The rise of the chamber's gas temperature since the start at nominal flux
`flux` [kW/m²], as a function of time [s] returning kelvins: the
calcium-silicate rear-face record at that flux (`silcal_data`, to 600 s),
continued past it at its mean slope over the last 200 s.
"""
function chamber_rise(flux; data = silcal_data())
    d = data[round(Int, flux)]
    T = d.T[:, end]
    ts, T0 = d.t, T[1]
    k = searchsortedfirst(ts, ts[end] - 200)
    slope = (T[end] - T[k]) / (ts[end] - ts[k])
    return t -> t <= ts[end] ? interp(ts, T, t) - T0 : T[end] - T0 + slope * (t - ts[end])
end

"""
    problem(species, flux, rho; T0, t_end, back = :chamber, T40 = nothing,
            h_front = nothing, cells = 100, save_every = 5.0, T_window = nothing,
            kwargs...)
        -> (problem, nodes)

The 1D block of 50 mm across the grain for `species` at nominal flux `flux`
[kW/m²] and dry density `rho` [kg/m³]. `T0` is the initial temperature [K],
a number or a `(depths, temperatures)` pair interpolated in depth.
- **Exposed face** (`:top`, z = L): absorbs the nominal flux at the
  material's surface absorptivity (the wood's, then char 0.95, by volume
  share), re-radiates to `T_window(t)` (default T_a(t)) and convects to T_a(t)
  at the thesis's h (chamber.toml [h_front]), or `h_front` when given.
- **Gas temperature** T_a(t) = T_ref + the calcium-silicate rear-face rise at
  this flux (`chamber_rise`), with T_ref the initial temperature at 40 mm.
- **Back face** (`:bottom`), by `back`:
  - `:chamber`, T = T_a(t), the thesis's condition;
  - `:tc40`, the domain ends at 40 mm, held at the measured 40 mm record `T40`
    (a `(t, T)` pair; temperature runs only);
  - `:adiabatic`, no flux.
- **Gas**: out through the exposed face to argon at T_a(t); the back is
  impermeable.
`kwargs` go to [`material`](@ref).
"""
function problem(
    species::AbstractString,
    flux::Real,
    rho::Real;
    T0,
    t_end::Real,
    back::Symbol = :chamber,
    T40 = nothing,
    cells::Int = 100,
    save_every::Real = 5.0,
    T_window = nothing,
    h_front::Union{Nothing,Real} = nothing,
    ch = chamber(),
    kwargs...,
)
    d = deck(species)
    mat = material(d, rho; ch, kwargs...)
    L = back === :tc40 ? 0.04 : Float64(ch["thickness"])
    T_init = T0 isa Real ? (z -> Float64(T0)) : (z -> interp(T0[1], T0[2], L - z))
    T_ref = T0 isa Real ? Float64(T0) : interp(T0[1], T0[2], 0.04)
    rise = chamber_rise(flux)
    T_a = t -> T_ref + rise(t)
    T_rad = something(T_window, T_a)
    q = 1000 * Float64(flux)
    h = Float64(something(h_front, ch["h_front"][string(round(Int, flux))]))
    top = RadiativeFluxBC(q) + RadiativeBC(T_inf = T_rad) + ConvectiveBC(h = h, T_inf = T_a)
    bottom = if back === :chamber
        DirichletThermalBC(T_a)
    elseif back === :tc40
        T40 === nothing && error("back = :tc40 needs the 40 mm record T40")
        DirichletThermalBC(t -> interp(T40[1], T40[2], t))
    elseif back === :adiabatic
        HeatFluxBC(0.0)
    else
        error("unknown back condition $back")
    end
    reservoir = GasReservoir(P = P_ATM, T = T_a, mole_fractions = (ARGON = 1.0,))
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = bottom),
        mass = (top = AmbientGasBC(reservoir), bottom = ImpermeableBC()),
    )
    α = d["kinetics"]["alpha"]
    ξ0 = (α[1] * rho, α[2] * rho, α[3] * rho, 0.0, 0.0, 0.0, 0.0, 0.0)
    ξinit = [
        (z -> ξ0[j] for j = 1:7)...,
        z -> background_gas_concentration(mat, :ARGON, ξ0, T_init(z)),
    ]
    nodes = collect(range(0.0, L, length = cells + 1))
    prob = PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = mat,
        bc_set = bcs,
        T_initial = T_init,
        ξ_initial = ξinit,
        tspan = (0.0, Float64(t_end)),
        dt_initial = 1e-3,
        output_times = collect(0.0:save_every:Float64(t_end)),
    )
    return prob, nodes
end

# The initial profile of a temperature run: its live thermocouples' first
# readings against depth [m].
function _initial_profile(r::TemperatureRun)
    j = [j for j in eachindex(DEPTHS) if live(r, j) && isfinite(r.T0[j])]
    return (DEPTHS[j], r.T0[j])
end

_T0(r::MassRun) = silcal_data()[round(Int, r.flux)].T[1, end]
_T0(r::TemperatureRun) = _initial_profile(r)
_t_end(r::MassRun) = r.t[end]
_t_end(r::TemperatureRun) = min(r.t_end, r.t[end])

"""
    problem(r::MassRun; kwargs...)
    problem(r::TemperatureRun; back = :chamber, kwargs...)

The test's problem: its species, flux and density, to the end of its record
or exposure. A mass test starts uniform at the chamber's initial gas
temperature at its flux, the first reading of the calcium-silicate rear face
(its own initial temperature was not recorded); a temperature test starts from
its first row.
"""
function problem(r::MassRun; kwargs...)
    return problem(r.species, r.flux, r.rho; T0 = _T0(r), t_end = _t_end(r), kwargs...)
end

function problem(r::TemperatureRun; back::Symbol = :chamber, kwargs...)
    T40 = nothing
    if back === :tc40
        j = length(DEPTHS)
        live(r, j) || error("$(r.name): no 40 mm record")
        T40 = (r.t, r.T[:, j])
    end
    return problem(
        r.species,
        r.flux,
        r.rho;
        T0 = _T0(r),
        t_end = _t_end(r),
        back,
        T40,
        kwargs...,
    )
end

"""
    simulate(r; reltol = 1e-6, dtmax = 5.0, kwargs...) -> NamedTuple

Solve the test's problem with a `ConservationRecorder` and return the series
the comparison reads: temperatures [K] at the nominal thermocouple depths
(cell centres interpolated linearly), the surface and back temperatures [K],
the condensed mass over its initial value, the mass loss rate [g/(m² s)],
and the mass and energy closures of `conservation_audit`.
"""
function simulate(r; reltol = 1e-6, dtmax = 5.0, kwargs...)
    prob, nodes = problem(r; kwargs...)
    rec = ConservationRecorder()
    sol = solve(prob; reltol, abstol = 1e-6, dtmax, callback = rec)
    audit = conservation_audit(sol, rec)
    L = nodes[end]
    nt = length(sol.t)
    centres = 0.5 .* (nodes[1:(end-1)] .+ nodes[2:end])
    depths = filter(d -> d < L, DEPTHS)
    T = [interp(centres, view(sol.T, :, k), L - d) for k = 1:nt, d in depths]
    solids = 1:6
    m = [
        sum(
            sol.ξ[j, c, k] * (nodes[c+1] - nodes[c]) for
            j in solids, c = 1:(length(nodes)-1)
        ) for k = 1:nt
    ]
    return (;
        t = sol.t,
        depths,
        T,
        T_surface = sol.surface_T,
        T_back = sol.extras.back_temperature,
        m = m ./ m[1],
        mlr = 1000 .* sol.mass_loss_rate,
        retcode = string(sol.retcode),
        solve_time = sol.solve_time,
        audit,
    )
end

end # module
