# Fredlund's wet-wood cubes (Lund, TVBB-1003, 1988; TVBB-3020, 1985): spruce
# at 14.5 % and pine at 24 % moisture, 135 mm cubes heated across the grain at
# 90 kW/m² nominal in air, with temperatures and gauge pore pressures at
# 10-70 mm; and a 150 mm pine cube at 14 % with resistive moisture sensors.
# The data are read from the reports' figures (data/, see README.md). The
# deck (deck.toml) is the input of Fredlund's WOOD1 model; the moisture model
# (the phase change of the bound water and its vapour) and the one-way law
# compared with it are those of scripts/lemta_moist_wood.
module FredlundWetWood

using Pyrolysis
using Printf
using Statistics
using TOML

include(joinpath(@__DIR__, "..", "lemta_moist_wood", "LemtaMoistWood.jl"))
using .LemtaMoistWood: LAWS, PHASE_CHANGE_LAWS, ISOTHERMS, D_VAPOUR_AIR, heat_of_desorption

export CASES,
    DRY,
    deck,
    temperatures,
    pressures,
    moisture_1985,
    material,
    problem,
    simulate,
    write_series,
    read_series,
    series_path,
    RESULTS,
    DEPTHS,
    interp,
    crossing

const HERE = @__DIR__
const DATA = joinpath(HERE, "data")
const RESULTS = joinpath(HERE, "results")
const P_ATM = 101325.0
const σ_SB = 5.670374419e-8
const MWG = 9.79e3   # Pa per metre of water gauge (F88)

"""The sensor depths below the exposed face [mm]."""
const DEPTHS = (10, 20, 30, 40, 50, 60, 70)

deck() = TOML.parsefile(joinpath(HERE, "deck.toml"))
const CASES = ("spruce14", "pine24", "pine14_1985")
const DRY = "spruce_dry"

"""Piecewise linear in `xs`, held constant outside them; NaN at NaN (a trial
state the solver rejects)."""
function interp(xs::AbstractVector, ys::AbstractVector, x::Real)
    isnan(x) && return x
    x <= xs[1] && return float(ys[1])
    x >= xs[end] && return float(ys[end])
    i = searchsortedlast(xs, x)
    w = (x - xs[i]) / (xs[i+1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i+1]
end

"""First time a (t, y) series reaches `level` (NaN if never)."""
function crossing(t, y, level)
    k = findfirst(>=(level), y)
    k === nothing && return NaN
    k == 1 && return t[1]
    return t[k-1] + (level - y[k-1]) / (y[k] - y[k-1]) * (t[k] - t[k-1])
end

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

_rows(name) = [
    split(l, ',') for l in
    filter(l -> !isempty(l) && !startswith(l, "#"), readlines(joinpath(DATA, name)))[2:end]
]
_num(s) = (s == "" || s == "NaN") ? NaN : parse(Float64, s)

"""Measured temperatures of `case` on the exposure clock: Dict(depth_mm =>
(t [s], T [K]))."""
function temperatures(case::AbstractString; d = deck())
    off = 60 * d["cases"][case]["t_offset_min"]
    out = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    for r in _rows("temperatures.csv")
        r[1] == case || continue
        t, T = get!(out, parse(Int, r[2]), (Float64[], Float64[]))
        push!(t, 60 * _num(r[3]) - off)
        push!(T, _num(r[4]) + 273.15)
    end
    return out
end

"""Measured peak gauge pressures of `case`: rows (depth_mm, p_peak [Pa],
t_from, t_to [s, exposure clock], status)."""
function pressures(case::AbstractString; d = deck())
    off = 60 * d["cases"][case]["t_offset_min"]
    return [
        (
            depth = parse(Int, r[2]),
            p = MWG * _num(r[3]),
            t_from = 60 * _num(r[4]) - off,
            t_to = 60 * _num(r[5]) - off,
            status = String(r[6]),
        ) for r in _rows("pressures.csv") if r[1] == case
    ]
end

"""The 1985 moisture record (F85 Fig. 5.4 against Fig. 5.3), one row per depth."""
moisture_1985() = [
    (
        depth = parse(Int, r[1]),
        X_peak = _num(r[2]) / 100,
        t_peak = 60 * _num(r[3]),
        T_at_peak = _num(r[4]) + 273.15,
        t_100C = 60 * _num(r[5]),
        t_back = 60 * _num(r[6]),
        t_zero = 60 * _num(r[7]),
    ) for r in _rows("moisture_1985.csv")
]

"""The measured 10 mm record of `case` on the exposure clock: (t [s], T [K])."""
function tc10(case::AbstractString; d = deck())
    off = 60 * d["cases"][case]["t_offset_min"]
    rows = filter(r -> r[1] == case, _rows("tc10.csv"))
    t = [max(60 * _num(r[2]) - off, 0.0) for r in rows]
    T = [_num(r[3]) + 273.15 for r in rows]
    keep = [i == length(t) || t[i] < t[i+1] for i in eachindex(t)]
    return t[keep], T[keep]
end

"""The incident flux of `history` (data/flux.csv) as a function of time [s] → W/m²."""
function flux(history::AbstractString)
    rows = filter(r -> r[1] == history, _rows("flux.csv"))
    t = [60 * _num(r[2]) for r in rows]
    q = [1e3 * _num(r[3]) for r in rows]
    return s -> interp(t, q, s)
end

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

_table(Ts, ys) = (x = Float64.(Ts) .+ 273.15, y = Float64.(ys))
_tablelaw(tb) = T -> interp(tb.x, tb.y, T)

"""Component indices of a material built by [`material`](@ref)."""
const IX = (wood = 1, char = 2, water = 3, vapour = 4, volatiles = 5, air = 6)

"""
    material(d, case; law = :sorption, λ_scale, κ_scale = 1, k_scale = 1,
             kc_scale = 1, exchange_rate = 1e3, tortuosity = nothing)

WOOD1's wood and char of `case` (ρ_dry, the char at 120 kg/m³ bulk holding
the wood's volume), bound water, its vapour, the volatiles and air. The wood
pyrolyses to char and volatiles by Fredlund's in-situ kinetics, gated off
below 417 K. Under `:sorption` or `:freewater` the water and its vapour
exchange two ways toward equilibrium on the deck's isotherm or as free water
at `exchange_rate` [1/s], the wood the sorbing solid and the char not; under
a one-way law of `LAWS` the water desorbs by that law with the LEMTA decks'
heat of desorption. `λ_scale` multiplies the wood's transfer coefficient (by
default the deck's `transfer_wet_scale` at `wet_moisture` or more),
`κ_scale` both solids' permeability, `k_scale` their conductivity and
`kc_scale` the char's alone; a `tortuosity` replaces the lumped
transfer by the dusty-gas closure. The water holds its own volume at 1000
kg/m³ and adds to the conductivity by that share at 0.6 W/(m K).
"""
function material(
    d,
    case;
    law::Symbol = :sorption,
    λ_scale::Union{Nothing,Real} = nothing,
    κ_scale::Real = 1.0,
    k_scale::Real = 1.0,
    kc_scale::Real = 1.0,
    exchange_rate::Real = 1e3,
    tortuosity::Union{Nothing,Real} = nothing,
)
    c = d["cases"][case]
    ρ = Float64(c["rho_dry"])
    X0 = Float64(c["moisture"])
    λs = something(λ_scale, X0 >= d["wet_moisture"] ? d["transfer_wet_scale"] : 1.0)
    γ = d["rho_char"] / ρ
    ρs = Float64(d["skeletal_density"])
    wood = SolidComponent(
        :WOOD;
        ρ = ρ,
        c = CallableProperty(_tablelaw(_table(d["c_wood_T"], d["c_wood"]))),
        k = CallableProperty(_tablelaw(_table(d["k_wood_T"], k_scale .* d["k_wood"]))),
        ε = d["wood_emissivity"],
        λ = λs * d["transfer_wood"],
        κ = κ_scale * d["kappa_wood"],
        ρ_intrinsic = ρs,
    )
    char = SolidComponent(
        :CHAR;
        ρ = Float64(d["rho_char"]),
        c = CallableProperty(_tablelaw(_table(d["c_char_T"], d["c_char"]))),
        k = CallableProperty(
            _tablelaw(_table(d["k_char_T"], (k_scale * kc_scale) .* d["k_char"])),
        ),
        ε = d["char_emissivity"],
        λ = d["transfer_char"],
        κ = κ_scale * d["kappa_char"],
        ρ_intrinsic = ρs,
    )
    water = LiquidComponent(
        :MOISTURE;
        ρ = 1000.0,
        c = d["c_water"],
        k = 0.6,
        ε = d["wood_emissivity"],
    )
    vapour = GaseousComponent(
        :WATER_VAPOR;
        M = 0.018015,
        c = CallableProperty(_tablelaw(_table(d["c_vapour_T"], d["c_vapour"]))),
        k = 0.025,
    )
    volatiles = GaseousComponent(
        :VOLATILES;
        M = d["volatiles_M"],
        c = d["volatiles_c"],
        k = d["volatiles_k"],
    )
    air = GaseousComponent(:AIR; M = d["air_M"], c = d["air_c"], k = d["air_k"])
    pyrolysis = Reaction(
        :pyrolysis,
        :WOOD => (:CHAR, :VOLATILES);
        A = d["A_pyrolysis"],
        E = d["E_pyrolysis"],
        h = 0.0,
        yields = (γ, 1 - γ),
        T_min = d["T_min_pyrolysis"],
        gate_width = 2.0,
    )
    evaporation, phase_changes = if law in PHASE_CHANGE_LAWS
        (),
        (
            PhaseChange(
                :moisture,
                :MOISTURE => :WATER_VAPOR;
                activity = law === :freewater ? FreeLiquid() : ISOTHERMS[d["isotherm"]],
                sorbents = (:WOOD,),
                exchange_rate,
            ),
        )
    else
        p = LAWS[law]
        (
            (
                Reaction(
                    :evaporation,
                    :MOISTURE => :WATER_VAPOR;
                    A = p.A,
                    E = p.E,
                    h = StateDependentProperty(
                        (T, ξ) -> heat_of_desorption(
                            T,
                            ξ[IX.water] / max(ξ[IX.wood] + ξ[IX.char], 1e-12),
                        ),
                    ),
                    n = 1.0,
                    T_min = p.T_min,
                    gate_width = p.width,
                ),
            ),
            (),
        )
    end
    return Material(
        name = Symbol(case),
        components = (wood, char, water, vapour, volatiles, air),
        reactions = (pyrolysis, evaporation...),
        phase_changes = phase_changes,
        mixing = MaterialMixing(
            k = MixingSpec(PARALLEL),
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

# ---------------------------------------------------------------------------
# Scenario
# ---------------------------------------------------------------------------

"""
    problem(case; law = :sorption, front, cells, flux_scale = 1, h_front,
            openback = false, save_every = 15.0, kwargs...)
        -> (problem, nodes, top_depth)

The cube as a 1D column across the grain, by `front` (default the case's):
- `:tc10`: the column below the 10 mm thermocouple (`top_depth` 0.010 m),
  its top held at that thermocouple's measured record (`tc10`); no flux
  enters;
- `:flux`: the whole column; `:top` absorbs `flux_scale` times the incident
  flux of the case's history at the surface absorptivity (by volume share),
  re-radiates at ε to `T_amb` and convects at `h_front`.
Gas leaves through the top to air at 1 atm whose vapour is the equilibrium
vapour of the starting wood at `T_amb` (the humidity that leaves the
conditioned wood at rest). The back is adiabatic and impermeable (`openback`: open to
the same air). Uniform start at T0 with the moisture's water; under a phase
change the pore vapour starts in equilibrium and the air makes up 1 atm.
The cells are 0.5 mm unless `cells` says otherwise. `kwargs` go to
[`material`](@ref).
"""
function problem(
    case::AbstractString;
    law::Symbol = :sorption,
    front::Union{Nothing,Symbol} = nothing,
    cells::Union{Nothing,Int} = nothing,
    flux_scale::Real = 1.0,
    h_front::Union{Nothing,Real} = nothing,
    openback::Bool = false,
    save_every::Real = 15.0,
    d = deck(),
    kwargs...,
)
    c = d["cases"][case]
    mat = material(d, case; law, kwargs...)
    fr = something(front, Symbol(c["front"]))
    top_depth = fr === :tc10 ? 0.010 : 0.0
    L = Float64(c["thickness"]) - top_depth
    n = something(cells, round(Int, L / 5e-4))
    Ta = Float64(d["T_amb"])
    T0 = Float64(get(c, "T0", d["T0"]))
    ρ = Float64(c["rho_dry"])
    X0 = Float64(c["moisture"])
    ρw = ρ / (1 + X0 * ρ / 1000)   # the water holds its own volume
    ξ0 = [ρw, 0.0, X0 * ρw, 0.0, 0.0, 0.0]
    if isempty(mat.phase_changes) || X0 == 0
        ξ0[IX.air] = background_gas_concentration(mat, :AIR, Tuple(ξ0), T0)
        x_v = 0.0
    else
        p_eq = equilibrium_vapour_pressure(mat, Tuple(ξ0), T0)
        ξ0[IX.vapour] = equilibrium_vapour_concentration(mat, Tuple(ξ0), T0)
        ξ0[IX.air] =
            background_gas_concentration(mat, :AIR, Tuple(ξ0), T0; P = P_ATM - p_eq)
        x_v = equilibrium_vapour_pressure(mat, Tuple(ξ0), Ta) / P_ATM
    end
    top = if fr === :tc10
        tt, TT = tc10(case; d)
        DirichletThermalBC(t -> interp(tt, TT, t))
    elseif fr === :flux
        q = flux(c["flux"])
        f = Float64(flux_scale)
        h = Float64(something(h_front, d["h_front"]))
        RadiativeFluxBC(t -> f * q(t)) +
        RadiativeBC(T_inf = Ta) +
        ConvectiveBC(h = h, T_inf = Ta)
    else
        error("unknown front $fr")
    end
    reservoir = GasReservoir(
        P = P_ATM,
        T = Ta,
        mole_fractions = (AIR = 1.0 - x_v, WATER_VAPOR = x_v),
    )
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = AdiabaticBC()),
        mass = (
            top = AmbientGasBC(reservoir),
            bottom = openback ? AmbientGasBC(reservoir) : ImpermeableBC(),
        ),
    )
    nodes = collect(range(0.0, L, length = n + 1))
    tend = Float64(c["t_end"])
    prob = PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = mat,
        bc_set = bcs,
        T_initial = z -> T0,
        ξ_initial = [(z -> ξ0[j]) for j = 1:6],
        tspan = (0.0, tend),
        dt_initial = 1e-2,
        output_times = collect(0.0:save_every:tend),
    )
    return prob, nodes, top_depth
end

"""
    simulate(case; law = :sorption, reltol = 1e-6, abstol = 1e-6,
             dtmax = 10.0, kwargs...) -> NamedTuple

Solve the case and return the series the comparison reads, at the sensor
depths `DEPTHS` below the exposed face (cell centres interpolated): the
temperature [K] (`T10` …), the gauge pore pressure [Pa] (`P10` …) and the
moisture content of the wood (`X10` …); with the surface temperature, the
mass per unit area [kg/m²], the water left, the largest gauge pressure, the
largest moisture content, and the evaporation and condensation rates
[kg/(m² s)] (the phase change's gross exchange; under a one-way law its
rate), and the mass and energy closures of `conservation_audit`. Under a
phase change `phase_change` holds the largest ratio of the pore vapour
pressure to its equilibrium pressure where the wood holds water
(X > 10⁻³), and the accepted and rejected steps.
"""
function simulate(
    case::AbstractString;
    law::Symbol = :sorption,
    reltol = 1e-6,
    abstol = 1e-6,
    dtmax = 10.0,
    kwargs...,
)
    d = deck()
    prob, nodes, top_depth = problem(case; law, d, kwargs...)
    rec = ConservationRecorder()
    sol = solve(prob; reltol, abstol, dtmax, callback = rec)
    audit = conservation_audit(sol, rec)
    mat = prob.material
    nt = length(sol.t)
    L = nodes[end]
    w = diff(nodes)
    zc = 0.5 .* (nodes[1:(end-1)] .+ nodes[2:end])
    NC = size(sol.ξ, 1)
    two_way = !isempty(mat.phase_changes)
    gv = findfirst(==(:WATER_VAPOR), sol.gas.species)
    # a sensor at `depth` below the exposed face; at or above the column's
    # top (the 10 mm boundary) the reading is the top cell's
    at(v, depth) = interp(zc, v, L - (depth / 1000 - top_depth))
    Td = Dict(dd => zeros(nt) for dd in DEPTHS)
    Pd = Dict(dd => zeros(nt) for dd in DEPTHS)
    Xd = Dict(dd => zeros(nt) for dd in DEPTHS)
    water = zeros(nt)
    pmax = zeros(nt)
    Xmax = zeros(nt)
    ratio = zeros(nt)
    evap = zeros(nt)
    rates = zeros(length(mat.reactions))
    heats = similar(rates)
    P = zeros(length(w))
    X = zeros(length(w))
    for k = 1:nt
        for c in eachindex(w)
            ξ = ntuple(j -> sol.ξ[j, c, k], NC)
            T = sol.T[c, k]
            P[c] = Pyrolysis.Physics.pore_state(ξ, T, mat).P - P_ATM
            X[c] = ξ[IX.wood] > 1e-9 ? ξ[IX.water] / ξ[IX.wood] : 0.0
            water[k] += ξ[IX.water] * w[c]
            if two_way
                X[c] > 1e-3 && (
                    ratio[k] = max(
                        ratio[k],
                        sol.gas.partial_pressure[gv, c, k] /
                        equilibrium_vapour_pressure(mat, ξ, T),
                    )
                )
            else
                Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, mat)
                evap[k] += rates[2] * w[c]
            end
        end
        pmax[k] = maximum(P)
        Xmax[k] = maximum(X)
        for dd in DEPTHS
            Td[dd][k] = at(view(sol.T, :, k), dd)
            Pd[dd][k] = at(P, dd)
            Xd[dd][k] = at(X, dd)
        end
    end
    cond = zeros(nt)
    if two_way
        evap .= sol.extras.evaporation_rate
        cond .= sol.extras.condensation_rate
    end
    cols = Dict{Symbol,Vector{Float64}}()
    for dd in DEPTHS
        cols[Symbol("T$dd")] = Td[dd]
        cols[Symbol("P$dd")] = Pd[dd]
        cols[Symbol("X$dd")] = Xd[dd]
    end
    return (;
        t = sol.t,
        cols,
        T_surface = sol.surface_T,
        mass = sol.extras.total_mass,
        water,
        p_max = pmax,
        X_max = Xmax,
        evaporation = evap,
        condensation = cond,
        retcode = string(sol.retcode),
        solve_time = sol.solve_time,
        audit,
        phase_change = two_way ?
                       (;
            pv_over_peq = ratio,
            accepted = sol.statistics.accepted,
            rejected = sol.statistics.rejected,
        ) : nothing,
    )
end

# ---------------------------------------------------------------------------
# Series files
# ---------------------------------------------------------------------------

series_path(variant, name) = joinpath(RESULTS, variant, name * ".csv")

const DEPTH_COLUMNS = Tuple(Symbol(p, dd) for p in ("T", "P", "X") for dd in DEPTHS)
const COLUMNS = (:T_surface, :mass, :water, :p_max, :X_max, :evaporation, :condensation)

"""
    write_series(path, s)

The series of [`simulate`](@ref), one row per saved time: `t`, the depth
columns, `COLUMNS` and, under a phase change, `pv_over_peq`; then a line
`# retcode=… solve_time=… mass_closure=… energy_closure=…` and, under a
phase change, `# phase_change accepted=… rejected=…`.
"""
function write_series(path, s)
    mkpath(dirname(path))
    pc = s.phase_change
    names = (:t, DEPTH_COLUMNS..., COLUMNS..., (pc === nothing ? () : (:pv_over_peq,))...)
    get_(c, k) =
        c === :t ? s.t[k] :
        c in DEPTH_COLUMNS ? s.cols[c][k] :
        c === :pv_over_peq ? pc.pv_over_peq[k] : getfield(s, c)[k]
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for k in eachindex(s.t)
            println(io, join((@sprintf("%.6g", get_(c, k)) for c in names), ","))
        end
        a = s.audit
        println(
            io,
            "# retcode=$(s.retcode) solve_time=$(round(s.solve_time, digits = 2)) " *
            "mass_closure=$(@sprintf("%.2e", a.mass.relative)) " *
            "energy_closure=$(@sprintf("%.2e", a.energy.relative))",
        )
        pc === nothing ||
            println(io, "# phase_change accepted=$(pc.accepted) rejected=$(pc.rejected)")
    end
    return path
end

"""
    read_series(path) -> NamedTuple

A series written by [`write_series`](@ref): its columns by name, `meta` (the
`# retcode` line) and `phase_change` (`nothing`, or the `# phase_change`
line's values).
"""
function read_series(path)
    lines = readlines(path)
    data = filter(l -> !startswith(l, "#"), lines)
    header = Symbol.(split(data[1], ','))
    rows = [parse.(Float64, split(l, ',')) for l in data[2:end]]
    cols = NamedTuple{Tuple(header)}(Tuple([r[j] for r in rows] for j in eachindex(header)))
    meta = filter(l -> startswith(l, "# retcode"), lines)
    pcline = filter(l -> startswith(l, "# phase_change"), lines)
    phase_change =
        isempty(pcline) ? nothing :
        NamedTuple(
            Symbol(a) => parse(Int, b) for
            (a, b) in (split(f, '=') for f in split(pcline[1])[3:end])
        )
    return merge(cols, (; meta = isempty(meta) ? "" : meta[1], phase_change))
end

end # module
