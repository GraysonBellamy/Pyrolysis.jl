# The Lee et al. 2025 pressure-treated pine case (Fire Mater. 49, 85-100):
# 9.0 mm disks at 6.2 % moisture in CAPA II under nitrogen at 29 and
# 63 kW/m², their bottom-face temperature read by an IR camera on a painted
# copper foil, and their mass loss rate. The data are digitised from the
# paper's Figs 6 and 7 (data/, see README.md). The deck (deck.toml) is the
# material of examples/08 with the boundary conditions of the paper's
# ThermaKin2Ds condition file; its water is the phase change of the bound
# water and its vapour, with the paper's one-way step as the comparison.
module LeePtwCapa

using Pyrolysis
using Printf
using Statistics
using TOML

export RUNS,
    Run,
    deck,
    bottom_temperature,
    mass_loss_rate,
    material,
    problem,
    simulate,
    write_series,
    read_series,
    series_path,
    RESULTS,
    LAWS,
    interp,
    crossing_seconds

const HERE = @__DIR__
const DATA = joinpath(HERE, "data")
const RESULTS = joinpath(HERE, "results")
const P_ATM = 101325.0
const σ_SB = 5.670374419e-8

"""The moisture laws: the phase change on the deck's isotherm (`:sorption`),
the same as free water (`:freewater`), and the paper's one-way step
(`:thermakin`)."""
const LAWS = (:sorption, :freewater, :thermakin)

deck() = TOML.parsefile(joinpath(HERE, "deck.toml"))

"""Piecewise linear in `xs`, held constant outside them."""
function interp(xs::AbstractVector, ys::AbstractVector, x::Real)
    x <= xs[1] && return float(ys[1])
    x >= xs[end] && return float(ys[end])
    i = searchsortedlast(xs, x)
    w = (x - xs[i]) / (xs[i+1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i+1]
end

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

function _csv(name)
    lines = filter(l -> !isempty(l) && !startswith(l, "#"), readlines(joinpath(DATA, name)))
    return [split(l, ',') for l in lines[2:end]]
end

_num(s) = isempty(s) ? NaN : parse(Float64, s)

"""
One set of tests: the set radiant flux [kW/m²] and the end of the plotted
record [s].
"""
struct Run
    name::String
    flux::Float64
    t_end::Float64
end

const RUNS = Dict("63" => Run("63", 63.0, 450.0), "29" => Run("29", 29.0, 1100.0))

"""
The measured bottom-face temperature of `run` (Fig. 6): (t [s], T [K],
band_lo, band_hi, model [K], scored::Bool), where `model` is the paper's
ThermaKin line and `scored` is false where the reading is flagged.
"""
function bottom_temperature(run::AbstractString)
    rows = filter(r -> r[1] == run, _csv("fig6_bottom_temperature.csv"))
    col(j) = [_num(r[j]) for r in rows]
    return (
        t = col(2),
        T = col(3),
        lo = col(4),
        hi = col(5),
        model = col(6),
        scored = [r[7] == "0" for r in rows],
    )
end

"""The measured mass loss rate of `run` (Fig. 7) [kg/(m² s)], as
[`bottom_temperature`](@ref) without the flag."""
function mass_loss_rate(run::AbstractString)
    rows = filter(r -> r[1] == run, _csv("fig7_mlr.csv"))
    col(j) = [_num(r[j]) for r in rows]
    return (t = col(2), MLR = col(3), lo = col(4), hi = col(5), model = col(6))
end

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

"""Component names and indices of a material built by [`material`](@ref)."""
function layout(d)
    solids = Tuple(Symbol.(d["solids"]))
    ns = length(solids)
    return (;
        solids,
        ns,
        copper = ns + 1,
        water = ns + 2,
        vapour = ns + 3,
        gases = (ns+4):(ns+7),
        n2 = ns + 8,
    )
end

"""
    material(d; law = :sorption, λ_wood, κ_scale = 1, exchange_rate = 1e3,
             tortuosity = nothing)

The treated pine of examples/08 (five solids, four pyrolysis steps, their
gases), the copper foil, bound water, its vapour and nitrogen. Under
`:sorption` or `:freewater` the water and its vapour exchange two ways toward
equilibrium on the deck's isotherm or as free water at `exchange_rate`
[1/s], the wood stages the sorbing solids and the char not; under
`:thermakin` the water leaves by the paper's one-way Arrhenius step with its
constant heat. `λ_wood` is the transfer coefficient of the stages that are
still wood (PTW, PTW1), `κ_scale` multiplies every solid's permeability, and
a `tortuosity` replaces the lumped transfer by the dusty-gas closure.
- The water holds its own volume at 1000 kg/m³ and adds to the conductivity
  by that share at 0.6 W/(m K) (examples/08).
- The copper carries the wood's transfer coefficient and permeability, so it
  adds no resistance of its own to the gas below the wood; the face beneath
  it is impermeable.
"""
function material(
    d;
    law::Symbol = :sorption,
    λ_wood::Real = d["transfer_wood"],
    κ_scale::Real = 1.0,
    exchange_rate::Real = 1e3,
    tortuosity::Union{Nothing,Real} = nothing,
)
    law in LAWS || error("unknown law $law")
    ix = layout(d)
    wood_like = (1, 2)   # PTW, PTW1
    solids = Tuple(
        SolidComponent(
            ix.solids[j];
            ρ = d["rho"][j],
            c = let c0 = d["c0"][j], c1 = d["c1"][j]
                CallableProperty(T -> c0 + c1 * T)
            end,
            k = d["k"][j],
            ε = d["emissivity"][j],
            λ = j in wood_like ? Float64(λ_wood) : d["transfer_char"],
            κ = κ_scale * (j == 5 ? d["kappa_char"] : d["kappa_wood"]),
            ρ_intrinsic = d["skeletal"][j],
        ) for j = 1:ix.ns
    )
    copper = SolidComponent(
        :COPPER;
        ρ = d["foil_density"],
        c = d["foil_c"],
        k = d["foil_k"],
        ε = d["foil_emissivity"],
        λ = Float64(λ_wood),
        κ = κ_scale * d["kappa_wood"],
        ρ_intrinsic = d["foil_density"],
    )
    water = LiquidComponent(
        :WATER;
        ρ = 1000.0,
        c = CallableProperty(T -> 5.2e3 - 6.7 * T + 1.1e-2 * T^2),
        k = 0.6,
        ε = d["emissivity"][1],
    )
    vapour = GaseousComponent(
        :WATER_VAPOR;
        M = 0.018015,
        c = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2),
        k = 0.025,
    )
    gases = Tuple(
        GaseousComponent(
            Symbol("PTW_gas$j");
            M = d["gas_M"][j],
            c = d["gas_c"][j],
            k = 0.020,
        ) for j = 1:4
    )
    n2 = GaseousComponent(:N2; M = d["n2_M"], c = d["n2_c"], k = d["n2_k"])
    pyrolysis = Tuple(
        Reaction(
            Symbol("ptw_step$j"),
            ix.solids[j] => (ix.solids[j+1], Symbol("PTW_gas$j"));
            A = d["A"][j],
            E = d["E"][j],
            h = d["h"][j],
            yields = (d["yield"][j], 1 - d["yield"][j]),
        ) for j = 1:4
    )
    evaporation, phase_changes = if law === :thermakin
        (
            (
                Reaction(
                    :water_evaporation,
                    :WATER => :WATER_VAPOR;
                    A = d["water_A"],
                    E = d["water_E"],
                    h = d["water_h"],
                ),
            ),
            (),
        )
    else
        (
            (),
            (
                PhaseChange(
                    :moisture,
                    :WATER => :WATER_VAPOR;
                    activity = law === :freewater ? FreeLiquid() :
                               getfield(Pyrolysis, Symbol(d["isotherm"])),
                    sorbents = Tuple(Symbol.(d["sorbents"])),
                    exchange_rate,
                ),
            ),
        )
    end
    w = Float64(d["mixing_weight"])
    return Material(
        name = :pressure_treated_pine,
        components = (solids..., copper, water, vapour, gases..., n2),
        reactions = (evaporation..., pyrolysis...),
        phase_changes = phase_changes,
        mixing = MaterialMixing(
            k = MixingSpec(WEIGHTED, w),
            λ = MixingSpec(WEIGHTED, w),
            κ = MixingSpec(WEIGHTED, w),
        ),
        gas_transport = tortuosity === nothing ? LumpedTransfer() :
                        DustyGas(;
            tortuosity,
            D_ref = 2.6e-5,   # water vapour in air at 298 K (m²/s)
            T_ref = 298.15,
            P_ref = P_ATM,
        ),
    )
end

"""
The pore gas over the condensed concentrations `ξ0` at `T`, (vapour, N₂)
[kg/m³ bulk]: under a phase change the vapour in equilibrium with the water
and the N₂ making up the ambient pressure; otherwise N₂ alone.
"""
function start_pore_gas(mat, ξ0, T)
    isempty(mat.phase_changes) &&
        return (0.0, background_gas_concentration(mat, :N2, ξ0, T))
    p_eq = equilibrium_vapour_pressure(mat, ξ0, T)
    return (
        equilibrium_vapour_concentration(mat, ξ0, T),
        background_gas_concentration(mat, :N2, ξ0, T; P = P_ATM - p_eq),
    )
end

# ---------------------------------------------------------------------------
# Scenario
# ---------------------------------------------------------------------------

"""The nitrogen over the exposed face [K] at time `t` [s] (deck `n2_history`
above `T_background`)."""
function n2_temperature(d)
    a1, b1, a2, b2 = d["n2_history"]
    Tb = Float64(d["T_background"])
    return t -> Tb + a1 * exp(b1 * t) + a2 * exp(b2 * t)
end

"""
The back face's surroundings for `run`: (T_rad(t), h(t), T_conv(t)), the
temperature the foil exchanges radiation with, the convection coefficient
and the gas temperature of the convection.
"""
function back_surroundings(r::Run, d)
    if r.name == "63"
        T1, T2, τ = d["back_T_start"], d["back_T_end"], d["back_ramp"]
        extra = σ_SB * (T2^4 - T1^4)
        ramp(t) = clamp(t / τ, 0.0, 1.0)
        return (
            t -> (T1^4 + extra * ramp(t) / σ_SB)^0.25,
            t -> d["h_back"] * ramp(t),
            t -> T2,
        )
    else
        ts, Ts = Float64.(d["back29_t"]), Float64.(d["back29_T"])
        wall(t) = interp(ts, Ts, t)
        return (wall, t -> Float64(d["h_back"]), wall)
    end
end

"""
    problem(r::Run; law = :sorption, cells = 180, n2 = :deck, foil = true,
            save_every = 2.0, kwargs...) -> (problem, nodes, ix)

The disk as a 1D column across its 9.0 mm and the copper foil below it:
- `:top` absorbs the set flux of `r` times the condition file's radial factor
  at the surface absorptivity (by volume share), gains σT_b⁴ from the
  background and re-radiates, and convects to the nitrogen (`n2 = :deck`, the
  condition file's history; `:background`, at T_b). Gas leaves to dry N₂ at
  1 atm (the reservoir names the vapour at zero).
- `:bottom` is the painted foil (ε 0.92), impermeable, exchanging radiation
  and convection with the sidewall (`back_surroundings`).
- `foil = true` puts the copper's mass in the bottom cells, at most half of
  each cell's volume, the wood filling the rest; the column is then 9.05 mm.
- Uniform start at the deck's T0 with the moisture's water; under a phase
  change the pore vapour starts in equilibrium and the N₂ makes up 1 atm.
`kwargs` go to [`material`](@ref).
"""
function problem(
    r::Run;
    law::Symbol = :sorption,
    cells::Int = 180,
    n2::Symbol = :deck,
    foil::Bool = true,
    save_every::Real = 2.0,
    t_end::Real = r.t_end,
    d = deck(),
    kwargs...,
)
    mat = material(d; law, kwargs...)
    ix = layout(d)
    Lw = Float64(d["thickness"])
    δ = Float64(d["foil_thickness"])
    L = foil ? Lw + δ : Lw
    Δ = L / cells
    nfoil = foil ? ceil(Int, 2δ / Δ) : 0
    f = foil ? δ / (nfoil * Δ) : 0.0   # copper's volume share in a foil cell
    z_foil = nfoil * Δ
    MC = Float64(d["moisture"])
    ρ0 = d["rho"][1]
    ρ_ptw = ρ0 / (1 + MC * ρ0 / 1000)   # the water holds its own volume
    nc = length(mat.components)
    function condensed(share_wood, share_cu)
        ξ = zeros(nc)
        ξ[1] = share_wood * ρ_ptw
        ξ[ix.copper] = share_cu * d["foil_density"]
        ξ[ix.water] = MC * ξ[1]
        return ξ
    end
    T0 = Float64(d["T0"])
    ξw = condensed(1.0, 0.0)
    ξf = condensed(1 - f, f)
    ξw[ix.vapour], ξw[ix.n2] = start_pore_gas(mat, Tuple(ξw), T0)
    ξf[ix.vapour], ξf[ix.n2] = start_pore_gas(mat, Tuple(ξf), T0)
    ξinit = Function[(z -> (z < z_foil ? ξf[j] : ξw[j])) for j = 1:nc]
    Tn2 = n2 === :deck ? n2_temperature(d) : (t -> Float64(d["T_background"]))
    q = 1e3 * r.flux * d["radial_factor"]
    top =
        RadiativeFluxBC(q) +
        RadiativeBC(T_inf = Float64(d["T_background"])) +
        ConvectiveBC(h = Float64(d["h_front"]), T_inf = Tn2)
    Trad, hb, Tconv = back_surroundings(r, d)
    bottom =
        RadiativeBC(ε = Float64(d["foil_emissivity"]), T_inf = Trad) +
        ConvectiveBC(h = hb, T_inf = Tconv)
    reservoir =
        GasReservoir(P = P_ATM, T = Tn2, mole_fractions = (N2 = 1.0, WATER_VAPOR = 0.0))
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = bottom),
        mass = (top = AmbientGasBC(reservoir), bottom = ImpermeableBC()),
    )
    nodes = collect(range(0.0, L, length = cells + 1))
    prob = PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = mat,
        bc_set = bcs,
        T_initial = z -> T0,
        ξ_initial = ξinit,
        tspan = (0.0, Float64(t_end)),
        dt_initial = 1e-3,
        output_times = collect(0.0:save_every:Float64(t_end)),
    )
    return prob, nodes, ix
end

"""
    simulate(r::Run; law = :sorption, reltol = 1e-6, abstol = 1e-6,
             dtmax = 5.0, kwargs...) -> NamedTuple

Solve the run and return the series the comparison reads: the bottom-face
temperature [K], the exposed face's, the mass loss rate and the mass per unit
area [kg/(m² s), kg/m²], the water left [kg/m²], the largest pore pressure
[Pa], the moisture content of the bottom wood cell and the largest over the
column, and the evaporation and condensation rates [kg/(m² s)] (the phase
change's gross exchange; under the one-way step its rate), with the mass and
energy closures of `conservation_audit`. Under a phase change `phase_change`
also holds the largest ratio of the pore vapour pressure to its equilibrium
pressure where the wood holds water (X > 10⁻³), and the accepted and
rejected steps.
"""
function simulate(
    r::Run;
    law::Symbol = :sorption,
    reltol = 1e-6,
    abstol = 1e-6,
    dtmax = 5.0,
    kwargs...,
)
    d = deck()
    prob, nodes, ix = problem(r; law, d, kwargs...)
    rec = ConservationRecorder()
    sol = solve(prob; reltol, abstol, dtmax, callback = rec)
    audit = conservation_audit(sol, rec)
    mat = prob.material
    nt = length(sol.t)
    w = diff(nodes)
    two_way = !isempty(mat.phase_changes)
    water = zeros(nt)
    pmax = fill(P_ATM, nt)
    Xmax = zeros(nt)
    Xback = zeros(nt)
    ratio = zeros(nt)
    evap = zeros(nt)
    rates = zeros(length(mat.reactions))
    heats = similar(rates)
    gv = findfirst(==(:WATER_VAPOR), sol.gas.species)
    NC = size(sol.ξ, 1)
    sorb = [findfirst(==(Symbol(s)), ix.solids) for s in d["sorbents"]]
    # the bottom cell that holds wood only (above the foil)
    cback = findfirst(c -> sol.ξ[ix.copper, c, 1] == 0, 1:size(sol.ξ, 2))
    for k = 1:nt
        for c in eachindex(w)
            ξ = ntuple(j -> sol.ξ[j, c, k], NC)
            T = sol.T[c, k]
            water[k] += ξ[ix.water] * w[c]
            pmax[k] = max(pmax[k], Pyrolysis.Physics.pore_state(ξ, T, mat).P)
            S = sum(max(ξ[j], 0.0) for j in sorb)
            X = S > 1e-9 ? ξ[ix.water] / S : 0.0
            Xmax[k] = max(Xmax[k], X)
            c == cback && (Xback[k] = X)
            if two_way
                X > 1e-3 && (
                    ratio[k] = max(
                        ratio[k],
                        sol.gas.partial_pressure[gv, c, k] /
                        equilibrium_vapour_pressure(mat, ξ, T),
                    )
                )
            else
                Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, mat)
                evap[k] += rates[1] * w[c]
            end
        end
    end
    cond = zeros(nt)
    if two_way
        evap .= sol.extras.evaporation_rate
        cond .= sol.extras.condensation_rate
    end
    phase_change =
        two_way ?
        (;
            pv_over_peq = ratio,
            accepted = sol.statistics.accepted,
            rejected = sol.statistics.rejected,
        ) : nothing
    return (;
        t = sol.t,
        T_bottom = sol.extras.back_temperature,
        T_surface = sol.surface_T,
        MLR = sol.extras.MLR_gauge,
        mass = sol.extras.total_mass,
        water,
        p_max = pmax,
        X_back = Xback,
        X_max = Xmax,
        evaporation = evap,
        condensation = cond,
        retcode = string(sol.retcode),
        solve_time = sol.solve_time,
        audit,
        phase_change,
    )
end

# ---------------------------------------------------------------------------
# Series files
# ---------------------------------------------------------------------------

series_path(variant, name) = joinpath(RESULTS, variant, name * ".csv")

const COLUMNS = (
    :t,
    :T_bottom,
    :T_surface,
    :MLR,
    :mass,
    :water,
    :p_max,
    :X_back,
    :X_max,
    :evaporation,
    :condensation,
)

const PHASE_CHANGE_COLUMNS = (:pv_over_peq,)

"""
    write_series(path, s)

The series of [`simulate`](@ref), one row per saved time (`COLUMNS`), and a
line `# retcode=… solve_time=… mass_closure=… energy_closure=…`. A phase
change adds the columns `PHASE_CHANGE_COLUMNS` and a line
`# phase_change accepted=… rejected=…`.
"""
function write_series(path, s)
    mkpath(dirname(path))
    pc = s.phase_change
    extra = pc === nothing ? () : PHASE_CHANGE_COLUMNS
    open(path, "w") do io
        println(io, join(string.((COLUMNS..., extra...)), ","))
        for k in eachindex(s.t)
            vals = (
                (getfield(s, c)[k] for c in COLUMNS)...,
                (getfield(pc, c)[k] for c in extra)...,
            )
            println(io, join((@sprintf("%.6g", v) for v in vals), ","))
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

"""First time [s] a (t, y) series reaches `level`; NaN if never."""
function crossing_seconds(t, y, level)
    k = findfirst(>=(level), y)
    k === nothing && return NaN
    k == 1 && return t[1]
    return t[k-1] + (level - y[k-1]) / (y[k] - y[k-1]) * (t[k] - t[k-1])
end

end # module
