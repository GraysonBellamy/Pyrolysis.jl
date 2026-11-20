# The Bilbao et al. 1996 wet-wood case (J. Anal. Appl. Pyrolysis 36, 81-97):
# Pinus pinaster slabs, 19 mm, at 0, 9 and 20 % moisture (dry basis), heated
# from one face by a radiant panel in N₂ under a time-varying flux, with
# temperatures at 1, 9 and 19 mm and the mass. The data are digitised from
# the paper's figures (data/, see README.md). The deck (deck.toml) is the
# paper's measured wood with the LEMTA decks' shared laws; the moisture model
# (the phase change of the bound water and its vapour) and the one-way
# evaporation laws compared with it are those of scripts/lemta_moist_wood.
module BilbaoWetWood

using Pyrolysis
using Printf
using Statistics
using TOML

include(joinpath(@__DIR__, "..", "lemta_moist_wood", "LemtaMoistWood.jl"))
using .LemtaMoistWood:
    LAWS, PHASE_CHANGE_LAWS, ISOTHERMS, D_VAPOUR_AIR, heat_of_desorption, p_sat
const c_water = LemtaMoistWood.c_water
const c_vapour = LemtaMoistWood.c_vapour

export RUNS,
    Run,
    deck,
    flux_table,
    flux,
    temperatures,
    conversion,
    material,
    problem,
    simulate,
    write_series,
    read_series,
    series_path,
    RESULTS,
    LAWS,
    PHASE_CHANGE_LAWS,
    interp,
    rms_by_depth,
    crossing_minutes,
    plateau_minutes,
    layout,
    SCHEMES

const HERE = @__DIR__
const DATA = joinpath(HERE, "data")
const RESULTS = joinpath(HERE, "results")
const P_ATM = 101325.0
const σ_SB = 5.670374419e-8

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
    header = split(lines[1], ',')
    return header, [split(l, ',') for l in lines[2:end]]
end

"""The incident flux of curves A–E [W/m²] on a 1-min grid: (t [s], Dict(curve => q))."""
function flux_table()
    header, rows = _csv("flux.csv")
    t = [60 * parse(Float64, r[1]) for r in rows]
    q = Dict(
        string(h[1]) => [1000 * parse(Float64, r[j]) for r in rows] for
        (j, h) in enumerate(header) if j > 1
    )
    return t, q
end

"""The incident flux of `curve` as a function of time [s] → W/m²."""
function flux(curve::AbstractString)
    t, q = flux_table()
    qc = q[curve]
    return s -> interp(t, qc, s)
end

"""
One test: its moisture content `Y0` (dry basis), flux curve, the end of its
record [s], and the figures it is read from.
"""
struct Run
    name::String
    Y0::Float64
    curve::String
    t_end::Float64
    figures::String
end

const RUNS = Dict(
    "dry_E" => Run("dry_E", 0.0, "E", 37 * 60.0, "Figs 3, 7"),
    "y09_E" => Run("y09_E", 0.09, "E", 40 * 60.0, "Figs 11, 7"),
    "y20_E" => Run("y20_E", 0.20, "E", 40 * 60.0, "Figs 12, 7"),
    "y20_E_repeat" => Run("y20_E_repeat", 0.20, "E", 66 * 60.0, "Fig. 4"),
    "y09_B" => Run("y09_B", 0.09, "B", 44 * 60.0, "Fig. 8"),
)

"""Measured temperatures of `run`: Dict(depth_mm => (t [s], T [K]))."""
function temperatures(run::AbstractString)
    _, rows = _csv("temperatures.csv")
    out = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    for r in rows
        r[1] == run || continue
        d = parse(Int, r[4])
        t, T = get!(out, d, (Float64[], Float64[]))
        push!(t, 60 * parse(Float64, r[3]))
        push!(T, parse(Float64, r[5]))
    end
    return out
end

"""Measured total conversion X_T of `run` (mass lost over initial mass): (t [s], X)."""
function conversion(run::AbstractString)
    _, rows = _csv("conversion.csv")
    t = Float64[]
    X = Float64[]
    for r in rows
        r[1] == run || continue
        push!(t, 60 * parse(Float64, r[3]))
        push!(X, parse(Float64, r[4]))
    end
    return t, X
end

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

_celsius_law(c) = T -> c[1] + c[2] * (T - 273.15)

k_wood(d, T; scale = 1.0) =
    scale *
    d["k_wood_scale"] *
    (d["k_wood"][1] + d["k_wood"][2] * (d["rho_dry"] / 1000) * sqrt(T / 433))

k_char(d, T) = (p = d["k_char"]; p["a"] * T + p["b"] + 16σ_SB * T^3 / (3p["beta"]))

"""
The pyrolysis schemes: `:bilbao`, one first-order step with the authors'
Arrhenius part and DSC heats (deck.toml); `:larch`, the three independent
pseudo-components fitted to larch TGA in N₂ (scripts/lemta_argon/decks/
larch.toml, Flity 2023 thesis Annex C), with the thesis's heats H1, H2 per g of
dry wood (H3 = 0; the heat-only char step C2 → C4 and its H4 are left out).
`alpha` the shares of the dry wood, `gamma` the char per kg of reactant,
`h` per kg of reactant [J/kg, positive endothermic].
"""
const SCHEMES = Dict(
    :larch => (
        A = [2.4e9, 1.3e12, 718.0],
        E = [127.7e3, 173.1e3, 84.9e3],
        alpha = [0.15, 0.61, 0.24],
        gamma = [0.0, 0.12, 0.60],
        h = [-46.1e3 / 0.15, 23.7e3 / 0.61, 0.0],
    ),
)

"""Component names and indices of a material built by [`material`](@ref)."""
function layout(kinetics::Symbol)
    woods = kinetics === :bilbao ? (:WOOD,) : (:WOOD_1, :WOOD_2, :WOOD_3)
    chars = kinetics === :bilbao ? (:CHAR,) : (:CHAR_2, :CHAR_3)
    nw, nc = length(woods), length(chars)
    n = nw + nc
    return (;
        woods,
        chars,
        nw,
        nc,
        solid = 1:n,
        water = n + 1,
        vapour = n + 2,
        volatiles = n + 3,
        n2 = n + 4,
    )
end

"""
Heat of the single-step pyrolysis per kg of wood consumed [J/kg, positive
endothermic]: the DSC heats per kg of mass lost (deck `heat_early` while the
local dry-basis conversion X_s ≤ `X_switch`, `heat_late` after) times the
volatile yield 1 − γ, switched smoothly over ±0.01 in X_s, with
X_s = (1 − γ)(1 − ξ_wood/ρ_dry).
"""
function pyrolysis_heat(d)
    γ = d["char_yield"]
    ρ = d["rho_dry"]
    h1, h2, Xs = d["heat_early"], d["heat_late"], d["X_switch"]
    return (T, ξ) -> begin
        X = (1 - γ) * (1 - ξ[1] / ρ)
        w = 0.5 * (1 + tanh((X - Xs) / 0.01))
        (1 - γ) * ((1 - w) * h1 + w * h2)
    end
end

"""Final char yield of a scheme [kg per kg of dry wood]."""
char_yield(d, kinetics) =
    kinetics === :bilbao ? Float64(d["char_yield"]) :
    sum(SCHEMES[kinetics].alpha .* SCHEMES[kinetics].gamma)

"""
    material(d; law, kinetics = :bilbao, shrinkage = nothing, κ_scale = 1,
             λ_scale = 1, sorption = 1, k_rule = PARALLEL, k_scale = 1,
             A_scale = 1, pyrolysis_scale = 1, exchange_rate = 1e3,
             tortuosity = nothing)

Dry wood (one component, or three pseudo-components), char, bound water,
water vapour, volatiles and N₂. The wood pyrolyses to char and volatiles by
`kinetics` (`SCHEMES`). Under a law of `PHASE_CHANGE_LAWS` the water and its
vapour exchange two ways toward equilibrium on the deck's isotherm (`:sorption`)
or as free water (`:freewater`) at the rate `exchange_rate` [1/s], the wood
components the sorbing solids and the char not; under a law of `LAWS` the
water desorbs one way by that rate law with the LEMTA decks' heat of
desorption (`sorption` scales its sorption part).
- Each wood component has the dry density, so they fill the volume in their
  shares. The char's component density makes a fully converted layer keep
  1 − `shrinkage` of the wood's volume (`nothing`: all of it, on a fixed
  mesh); a solve with shrinkage moves the mesh with the solid.
- The water holds its own volume at 1000 kg/m³ and, under `k_rule =
  PARALLEL`, adds to the conductivity by that share at 0.6 W/(m K) (the LEMTA
  convention); emissivity and absorptivity mix by volume share.
- `κ_scale` multiplies both solids' permeability, `λ_scale` the wood's
  transfer coefficient, `k_scale` the wood's and the char's conductivity,
  `A_scale` a one-way law's pre-exponential factor and `pyrolysis_scale` the
  pyrolysis ones. A `tortuosity` replaces the lumped transfer by the
  dusty-gas closure (as `LemtaMoistWood.material`).
"""
function material(
    d;
    law::Symbol,
    kinetics::Symbol = :bilbao,
    shrinkage::Union{Nothing,Real} = nothing,
    κ_scale::Real = 1.0,
    λ_scale::Real = 1.0,
    sorption::Real = 1.0,
    k_rule = PARALLEL,
    k_scale::Real = 1.0,
    A_scale::Real = 1.0,
    pyrolysis_scale::Real = 1.0,
    exchange_rate::Real = 1e3,
    tortuosity::Union{Nothing,Real} = nothing,
)
    ix = layout(kinetics)
    ρ = Float64(d["rho_dry"])
    γ = char_yield(d, kinetics)
    εw = Float64(d["wood_emissivity"])
    εc = Float64(d["char_emissivity"])
    c_w = CallableProperty(_celsius_law(d["c_wood"]))
    k_w = CallableProperty(T -> k_wood(d, T; scale = k_scale))
    wood(name) = SolidComponent(
        name;
        ρ = ρ,
        c = c_w,
        k = k_w,
        ε = εw,
        λ = λ_scale * d["transfer_wood"],
        κ = κ_scale * d["kappa_wood"],
        ρ_intrinsic = d["skeletal_density"],
    )
    # Each char holds its reactant's volume, times f: a fully converted layer
    # keeps 1 − shrinkage of its volume. Wood components that leave no char
    # (S1 of the three-component schemes) vanish, so the chars of the others
    # take the share 1 − shrinkage of the whole: f = (1 − s)/Σ(α of charring
    # components).
    s = shrinkage === nothing ? 0.0 : Float64(shrinkage)
    yields, shares = if kinetics === :bilbao
        [γ], [1.0]
    else
        SCHEMES[kinetics].gamma[2:3], SCHEMES[kinetics].alpha[2:3]
    end
    f = (1 - s) / sum(shares)
    c_c = CallableProperty(_celsius_law(d["c_char"]))
    k_c = CallableProperty(T -> k_scale * k_char(d, T))
    chars = [
        SolidComponent(
            name;
            ρ = y * ρ / f,
            c = c_c,
            k = k_c,
            ε = εc,
            λ = d["transfer_char"],
            κ = κ_scale * d["kappa_char"],
            ρ_intrinsic = d["skeletal_density"],
        ) for (name, y) in zip(ix.chars, yields)
    ]
    water = LiquidComponent(
        :MOISTURE;
        ρ = 1000.0,
        c = CallableProperty(c_water),
        k = 0.6,
        ε = εw,
    )
    vapour = GaseousComponent(
        :WATER_VAPOR;
        M = 0.018015,
        c = CallableProperty(c_vapour),
        k = 0.025,
    )
    volatiles = GaseousComponent(
        :VOLATILES;
        M = d["volatiles_M"],
        c = d["volatiles_c"],
        k = d["volatiles_k"],
    )
    n2 = GaseousComponent(:N2; M = d["n2_M"], c = d["n2_c"], k = d["n2_k"])
    iw, solid = ix.water, ix.solid
    evaporation, phase_changes = if law in PHASE_CHANGE_LAWS
        (),
        (
            PhaseChange(
                :moisture,
                :MOISTURE => :WATER_VAPOR;
                activity = law === :freewater ? FreeLiquid() : ISOTHERMS[d["isotherm"]],
                sorbents = ix.woods,
                exchange_rate,
            ),
        )
    else
        p = LAWS[law]
        (
            Reaction(
                :evaporation,
                :MOISTURE => :WATER_VAPOR;
                A = p.A * A_scale,
                E = p.E,
                h = StateDependentProperty(
                    (T, ξ) -> heat_of_desorption(
                        T,
                        ξ[iw] / max(sum(ξ[j] for j in solid), 1e-12);
                        sorption,
                    ),
                ),
                n = 1.0,
                T_min = p.T_min,
                gate_width = p.width,
            ),
        ),
        ()
    end
    pyrolysis = if kinetics === :bilbao
        (
            Reaction(
                :pyrolysis,
                :WOOD => (:CHAR, :VOLATILES);
                A = pyrolysis_scale * d["A_pyrolysis"],
                E = d["E_pyrolysis"],
                h = StateDependentProperty(pyrolysis_heat(d)),
                yields = (γ, 1 - γ),
            ),
        )
    else
        sch = SCHEMES[kinetics]
        Tuple(
            sch.gamma[j] == 0 ?
            Reaction(
                Symbol("R$j"),
                ix.woods[j] => :VOLATILES;
                A = pyrolysis_scale * sch.A[j],
                E = sch.E[j],
                h = sch.h[j],
            ) :
            Reaction(
                Symbol("R$j"),
                ix.woods[j] => (ix.chars[j-1], :VOLATILES);
                A = pyrolysis_scale * sch.A[j],
                E = sch.E[j],
                h = sch.h[j],
                yields = (sch.gamma[j], 1 - sch.gamma[j]),
            ) for j = 1:3
        )
    end
    return Material(
        name = :pinus_pinaster,
        components = (map(wood, ix.woods)..., chars..., water, vapour, volatiles, n2),
        reactions = (pyrolysis..., evaporation...),
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
The pore gas at the start, (vapour, N₂) [kg/m³ bulk], over the condensed
concentrations `ξ0` at `T`: under a phase change the vapour in equilibrium
with the water and the N₂ making up the ambient pressure, otherwise N₂ alone
at the ambient pressure.
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

"""
The run's measured 1 mm record as a function of time [s] → K, linear between
the readings, starting from the deck's T0 at t = 0.
"""
function tc1_record(r::Run; d = deck())
    meas = temperatures(r.name)
    haskey(meas, 1) || error("$(r.name) has no 1 mm record")
    t, T = meas[1]
    t = [0.0; t]
    T = [Float64(d["T0"]); T]
    return s -> interp(t, T, s)
end

"""
    problem(r::Run; front = :flux, law = :sorption, kinetics = :bilbao,
            flux_factor = 1.0, cells = 100, t_end, h_front, h_back = 0, T_amb,
            λ_scale, save_every = 10.0, kwargs...) -> (problem, nodes)

The 1D slab across the grain, by `front`:
- `:flux`: the whole 19 mm; the exposed face (`:top`) absorbs `flux_factor`
  times the incident flux of the run's curve at the surface absorptivity
  (wood 0.78, char 0.95, by volume share), re-radiates at ε to `T_amb` and,
  with `h_front` > 0, convects to it;
- `:tc1`: the 18 mm below the 1 mm thermocouple, its top held at the
  thermocouple's measured record (`tc1_record`); no flux enters the model;
- `:top1`: the 1 mm above the thermocouple, both faces at that record (the
  companion of `:tc1`, for the mass of the layer it leaves out).
The back (`:bottom`) is adiabatic, or with `h_back` > 0 loses heat to `T_amb`
through the fibreglass at that h. Gas leaves through the top to dry N₂ at
1 atm and `T_amb` (the reservoir names the vapour at zero); the back is
impermeable. Uniform start at the deck's T0 with Y0 ρ_dry of water; under a
phase change the pore vapour starts in equilibrium with it and the N₂ makes
up the ambient pressure. `λ_scale` multiplies the wood's transfer
coefficient, by default the deck's `transfer_wet_scale` for a run at
`wet_moisture` or more and one below. `kwargs` go to [`material`](@ref).
"""
function problem(
    r::Run;
    front::Symbol = :flux,
    law::Symbol = :sorption,
    kinetics::Symbol = :bilbao,
    flux_factor::Real = 1.0,
    cells::Int = 100,
    t_end::Real = r.t_end,
    h_front::Union{Nothing,Real} = nothing,
    h_back::Real = 0.0,
    T_amb::Union{Nothing,Real} = nothing,
    λ_scale::Union{Nothing,Real} = nothing,
    save_every::Real = 10.0,
    d = deck(),
    kwargs...,
)
    λs = something(λ_scale, r.Y0 >= d["wet_moisture"] ? d["transfer_wet_scale"] : 1.0)
    mat = material(d; law, kinetics, λ_scale = λs, kwargs...)
    ix = layout(kinetics)
    L = Float64(d["thickness"])
    Ta = Float64(something(T_amb, d["T_amb"]))
    h = Float64(something(h_front, d["h_front"]))
    bottom = h_back > 0 ? ConvectiveBC(h = Float64(h_back), T_inf = Ta) : HeatFluxBC(0.0)
    if front === :flux
        q = flux(r.curve)
        f = Float64(flux_factor)
        top = RadiativeFluxBC(t -> f * q(t)) + RadiativeBC(T_inf = Ta)
        h > 0 && (top = top + ConvectiveBC(h = h, T_inf = Ta))
    elseif front === :tc1
        L -= 0.001
        top = DirichletThermalBC(tc1_record(r; d))
    elseif front === :top1
        L = 0.001
        cells = 4
        top = DirichletThermalBC(tc1_record(r; d))
        bottom = DirichletThermalBC(tc1_record(r; d))
    else
        error("unknown front $front")
    end
    reservoir =
        GasReservoir(P = P_ATM, T = Ta, mole_fractions = (N2 = 1.0, WATER_VAPOR = 0.0))
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = bottom),
        mass = (top = AmbientGasBC(reservoir), bottom = ImpermeableBC()),
    )
    ρ = Float64(d["rho_dry"])
    T0 = Float64(d["T0"])
    shares = kinetics === :bilbao ? [1.0] : SCHEMES[kinetics].alpha
    ξ0 = (
        (shares .* ρ)...,
        ntuple(_ -> 0.0, ix.nc)...,   # chars
        r.Y0 * ρ,            # water
        0.0,                 # vapour, set below
        0.0,                 # volatiles
        0.0,                 # N₂, set below
    )
    ξ_v, ξ_n2 = start_pore_gas(mat, ξ0, T0)
    ξinit = Function[z -> ξ0[j] for j = 1:(ix.n2-1)]
    ξinit[ix.vapour] = z -> ξ_v
    push!(ξinit, z -> ξ_n2)
    nodes = collect(range(0.0, L, length = cells + 1))
    prob = PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = mat,
        bc_set = bcs,
        T_initial = z -> T0,
        ξ_initial = ξinit,
        tspan = (0.0, Float64(t_end)),
        dt_initial = 1e-2,
        output_times = collect(0.0:save_every:Float64(t_end)),
    )
    return prob, nodes
end

"""
    simulate(r::Run; law = :sorption, kinetics = :bilbao, shrinkage = nothing,
             reltol = 1e-6, abstol = 1e-6, dtmax = 10.0, kwargs...) -> NamedTuple

Solve the run and return the series the comparison reads: temperatures [K]
at the thermocouples' material points initially 1 and 9 mm deep (cell
centres interpolated) and at 19 mm (the back face), the surface temperature,
the slab thickness [m], the total conversion X_T = 1 − m/m₀ of the condensed
mass (wood, char, water), the water left [kg/m²], the largest pore pressure
[Pa] and vapour pressure over p_sat(T), the evaporation rate [kg/(m² s)] and
the share of it where the vapour exceeds saturation, with the mass and energy
closures of `conservation_audit`. With `shrinkage` the mesh moves with the
solid (`use_ale`).

Under a one-way law the evaporation is the reaction's rate. Under a phase
change it is the gross evaporation of the solution's extras, and
`phase_change` holds what the phase change adds: the gross condensation
[kg/(m² s)], the largest ratio of the pore vapour pressure to its equilibrium
pressure where the wood holds water (X > 10⁻³), the largest moisture
content, and the accepted and rejected steps.
"""
function simulate(
    r::Run;
    front::Symbol = :flux,
    law::Symbol = :sorption,
    kinetics::Symbol = :bilbao,
    shrinkage::Union{Nothing,Real} = nothing,
    reltol = 1e-6,
    abstol = 1e-6,
    dtmax = 10.0,
    kwargs...,
)
    prob, nodes = problem(r; front, law, kinetics, shrinkage, kwargs...)
    ale = shrinkage !== nothing
    ix = layout(kinetics)
    rec = ConservationRecorder()
    sol = solve(prob; reltol, abstol, dtmax, callback = rec, use_ale = ale)
    audit = conservation_audit(sol, rec)
    mat = prob.material
    L = nodes[end]
    nt = length(sol.t)
    z(k) = ale ? view(sol.z, :, k) : nodes
    centres(k) = 0.5 .* (z(k)[1:(end-1)] .+ z(k)[2:end])
    at(k, depth) = interp(centres(k), view(sol.T, :, k), interp(nodes, z(k), L - depth))
    top_depth = front === :tc1 ? 0.001 : 0.0   # the domain's top below the heated face
    T1 = front === :tc1 ? tc1_record(r).(sol.t) : [at(k, 0.001) for k = 1:nt]
    T9 = [at(k, 0.009 - top_depth) for k = 1:nt]
    T19 = sol.extras.back_temperature
    NC = size(sol.ξ, 1)
    solids = (ix.solid..., ix.water)
    m = zeros(nt)
    water = zeros(nt)
    pmax = fill(P_ATM, nt)
    smax = zeros(nt)
    evap = zeros(nt)
    evap_super = zeros(nt)
    rates = zeros(length(mat.reactions))
    heats = similar(rates)
    two_way = !isempty(mat.phase_changes)
    cond = zeros(nt)
    ratio = zeros(nt)
    Xmax = zeros(nt)
    gv = LemtaMoistWood.vapour_gas(sol)
    for k = 1:nt
        w = diff(z(k))
        for c in eachindex(w)
            ξ = ntuple(j -> sol.ξ[j, c, k], NC)
            m[k] += sum(ξ[j] for j in solids) * w[c]
            water[k] += ξ[ix.water] * w[c]
            T = sol.T[c, k]
            pore = Pyrolysis.Physics.pore_state(ξ, T, mat)
            pmax[k] = max(pmax[k], pore.P)
            s = pore.p[ix.vapour] / p_sat(T)
            smax[k] = max(smax[k], s)
            if two_way
                pc = only(mat.phase_changes)
                e = max(LemtaMoistWood.cell_rate(mat, sol, ξ, c, k), 0.0) * w[c]
                X = Pyrolysis.Physics.moisture_content(pc, ξ, map(x -> max(x, 0.0), ξ))
                Xmax[k] = max(Xmax[k], X)
                X > 1e-3 && (
                    ratio[k] = max(
                        ratio[k],
                        sol.gas.partial_pressure[gv, c, k] /
                        equilibrium_vapour_pressure(mat, ξ, T),
                    )
                )
            else
                Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, mat)
                e = rates[end] * w[c]
                evap[k] += e
            end
            s > 1 && (evap_super[k] += e)
        end
    end
    if two_way
        evap .= sol.extras.evaporation_rate
        cond .= sol.extras.condensation_rate
    end
    phase_change =
        two_way ?
        (;
            condensation = cond,
            pv_over_peq = ratio,
            X_max = Xmax,
            accepted = sol.statistics.accepted,
            rejected = sol.statistics.rejected,
        ) : nothing
    return (;
        t = sol.t,
        T1,
        T9,
        T19,
        T_surface = sol.surface_T,
        thickness = [z(k)[end] - z(k)[1] for k = 1:nt] .+ top_depth,
        X_T = front === :tc1 ?
              _with_top_layer(
            r,
            sol.t,
            m;
            law,
            kinetics,
            shrinkage,
            reltol,
            abstol,
            dtmax,
            kwargs...,
        ) : 1 .- m ./ m[1],
        water,
        p_max = pmax,
        supersaturation = smax,
        evaporation = evap,
        evaporation_supersaturated = evap_super,
        retcode = string(sol.retcode),
        solve_time = sol.solve_time,
        audit,
        phase_change,
    )
end

"""
X_T of a `front = :tc1` solve with the 1 mm layer above the thermocouple
added: that layer solved on its own at the thermocouple's record (`front =
:top1`, a lower bound on its temperature), its condensed mass added to the
domain's `m` [kg/m²] at the times `t`.
"""
function _with_top_layer(
    r,
    t,
    m;
    law,
    kinetics,
    shrinkage,
    reltol,
    abstol,
    dtmax,
    kwargs...,
)
    kw = Dict{Symbol,Any}(kwargs)
    pop!(kw, :cells, nothing)
    prob, _ = problem(r; front = :top1, law, kinetics, shrinkage, kw...)
    sol = solve(prob; reltol, abstol, dtmax, use_ale = shrinkage !== nothing)
    ix = layout(kinetics)
    solids = (ix.solid..., ix.water)
    zt(k) =
        shrinkage === nothing ? collect(range(0.0, 0.001, length = 5)) : view(sol.z, :, k)
    mt = [
        sum(sum(sol.ξ[j, c, k] for j in solids) * (zt(k)[c+1] - zt(k)[c]) for c = 1:4)
        for k in eachindex(sol.t)
    ]
    mtop = [interp(sol.t, mt, x) for x in t]
    total = m .+ mtop
    return 1 .- total ./ total[1]
end

# ---------------------------------------------------------------------------
# Series files
# ---------------------------------------------------------------------------

series_path(variant, name) = joinpath(RESULTS, variant, name * ".csv")

const COLUMNS = (
    :t,
    :T1,
    :T9,
    :T19,
    :T_surface,
    :thickness,
    :X_T,
    :water,
    :p_max,
    :supersaturation,
    :evaporation,
    :evaporation_supersaturated,
)

# The columns a phase change adds, from `s.phase_change`.
const PHASE_CHANGE_COLUMNS = (:condensation, :pv_over_peq, :X_max)

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
            println(
                io,
                join(
                    (
                        @sprintf("%.6g", v) for v in (
                            (getfield(s, c)[k] for c in COLUMNS)...,
                            (getfield(pc, c)[k] for c in extra)...,
                        )
                    ),
                    ",",
                ),
            )
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
`# retcode` line) and `phase_change`, `nothing` for a one-way law and
otherwise the `# phase_change` line's values (`accepted`, `rejected`).
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

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

"""
RMS of model minus measurement [K] at each depth (1, 9, 19 mm) over the
measured points within `window` [s] (model interpolated to the measured times).
"""
function rms_by_depth(run::AbstractString, s; window = (0.0, Inf))
    meas = temperatures(run)
    model = Dict(1 => s.T1, 9 => s.T9, 19 => s.T19)
    out = Dict{Int,Float64}()
    for (d, (t, T)) in meas
        keep = findall(x -> window[1] <= x <= min(window[2], s.t[end]), t)
        isempty(keep) && continue
        out[d] = sqrt(mean((interp(s.t, model[d], t[i]) - T[i])^2 for i in keep))
    end
    return out
end

"""First time [min] a (t [s], T [K]) series reaches `T_level`; NaN if never."""
function crossing_minutes(t, T, T_level)
    k = findfirst(>=(T_level), T)
    k === nothing && return NaN
    k == 1 && return t[1] / 60
    return (t[k-1] + (T_level - T[k-1]) / (T[k] - T[k-1]) * (t[k] - t[k-1])) / 60
end

"""Minutes a (t [s], T [K]) series spends within [lo, hi] before first leaving it upward."""
function plateau_minutes(t, T; lo = 365.0, hi = 400.0)
    a = crossing_minutes(t, T, lo)
    b = crossing_minutes(t, T, hi)
    return b - a
end

end # module
