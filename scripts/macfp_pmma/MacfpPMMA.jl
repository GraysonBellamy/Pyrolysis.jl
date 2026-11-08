# MaCFP PMMA decks for Pyrolysis.jl: canonical deck files (decks/*.toml) to
# problems, the MaCFP-2 target scenarios (1D gasification, TGA), and the
# MaCFP output quantities. See README.md.
module MacfpPMMA

using Pyrolysis
using TOML
using Printf

export apparatus_problem,
    apparatus_series,
    APPARATUS_TESTS,
    APPARATUS_HEATER,
    read_experiment,
    experiment_path,
    inert_problem,
    inert_series,
    INERT_TESTS,
    inert_experiment_path,
    load_deck,
    load_decks,
    material,
    gasification_problem,
    tga_problem,
    gasification_series,
    tga_series,
    write_series,
    read_series,
    read_matl_csv,
    matl_db,
    check_matl_db,
    GASIFICATION_CASES,
    TGA_RATES,
    prediction_path,
    tga_prediction_path,
    SOLVE_KWARGS,
    RESULTS,
    load_run,
    load_submission,
    metrics,
    top_column,
    burnout_time,
    crossing,
    interp,
    APPARATUS_PANELS,
    apparatus_experiment,
    group_predictions,
    nrms,
    apparatus_nrms,
    COST_TOL,
    REFLEVEL,
    combined_error,
    cost_rows,
    span_scaled

const HERE = @__DIR__
const DECKS = joinpath(HERE, "decks")
const T_AMBIENT = 293.0
const SAMPLE_AREA = 0.01          # m², 10 cm × 10 cm (MaCFP-2 output scaling)
const GAS_M = 0.10012             # kg/mol, MMA
const PHI = 0.01                  # nominal pore fraction of every solid
const LAMBDA = 2e-5               # lumped gas transfer coefficient of every solid
# The matl-db commit whose PMMA tree these decks and comparisons read.
const MATL_DB_COMMIT = "2772b0ebd81bfba727020550bbe0f5328aadccca"

# (incident flux kW/m², thickness mm, end time s): the MaCFP-2 target matrix.
# End times cover the submissions' spans (BUW-FZJ's 12000 s runs aside, whose
# inert residue never burns out).
const GASIFICATION_CASES = [
    (10, 6, 8000.0),
    (10, 12, 12000.0),
    (25, 6, 2000.0),
    (25, 12, 2400.0),
    (65, 6, 1500.0),
    (65, 12, 2000.0),
]
const TGA_RATES = (10, 100)         # K/min

matl_db() = get(ENV, "MATL_DB", joinpath(homedir(), "Documents", "git", "matl-db"))

"""Check that the matl-db checkout's PMMA tree is the pinned one."""
function check_matl_db()
    root = matl_db()
    isdir(joinpath(root, "PMMA")) || error("matl-db not found at $root (set MATL_DB)")
    pinned = readchomp(`git -C $root rev-parse $(MATL_DB_COMMIT):PMMA`)
    current = readchomp(`git -C $root rev-parse HEAD:PMMA`)
    pinned == current ||
        @warn "matl-db PMMA tree differs from the pinned commit" pinned current
    return root
end

load_deck(path::AbstractString) = TOML.parsefile(path)
function load_decks(; gasification_only = false)
    decks = [
        load_deck(joinpath(DECKS, f)) for f in sort(readdir(DECKS)) if endswith(f, ".toml")
    ]
    gasification_only && filter!(d -> !get(d, "kinetics_only", false), decks)
    return decks
end

# ---------------------------------------------------------------------------
# Property forms → functions of T [K]
# ---------------------------------------------------------------------------

function law(spec::Dict)
    form = spec["form"]
    if form == "constant"
        v = Float64(spec["value"])
        return T -> v
    elseif form == "linear"
        a, b = Float64(spec["a"]), Float64(spec["b"])
        return T -> a + b * T
    elseif form == "power"
        c0, Tr, n = Float64(spec["c0"]), Float64(spec["Tr"]), Float64(spec["n"])
        return T -> c0 * (T / Tr)^n
    elseif form == "table"
        Ts, vs = Float64.(spec["T"]), Float64.(spec["v"])
        return T -> _table(Ts, vs, T)
    elseif form == "switch"
        Tsw = Float64(spec["T"])
        lo, hi = law(spec["below"]), law(spec["above"])
        return T -> T < Tsw ? lo(T) : hi(T)
    else
        error("unknown property form $form")
    end
end

# Piecewise linear, held constant outside the table (FDS RAMP convention).
function _table(Ts, vs, T)
    T <= Ts[1] && return vs[1]
    T >= Ts[end] && return vs[end]
    i = searchsortedlast(Ts, T)
    w = (T - Ts[i]) / (Ts[i+1] - Ts[i])
    return (1 - w) * vs[i] + w * vs[i+1]
end

# Pyrolysis.jl requires k and c positive on [200, 3000] K. A law that turns
# non-positive inside that bracket is floored at 1 % of its value at 300 K:
# UMD's k = 0.34 − 4.2e-4 T above 395 K reaches zero at 810 K, which the
# surface-absorption variant passes late in its 65 kW/m² runs.
function positive(f, name)
    floor_value = 0.01 * abs(f(300.0))
    bad = any(T -> f(T) <= 0, 200.0:1.0:3000.0)
    bad && @info "property $name floored at $floor_value above its zero crossing"
    return bad ? (T -> max(f(T), floor_value)) : f
end

# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

"""
    material(deck; substrate = :none)

The deck's material: its condensed components, one MMA gas, its reactions.
An in-depth deck's components absorb and emit at their absorption coefficient
κ [1/m] (mass coefficient α = κ/ρ); a surface deck's absorb nothing in depth.
`substrate = :kaowool` adds the apparatus's Kaowool board as two components: an
opaque `:kaowool` for the cell under the PMMA and a transparent `:kaowool_bulk`
beneath it.
"""
function material(deck; substrate = :none)
    comps = deck["components"]
    gas_cp = Float64(get(deck["gas"], "cp", 0.0))
    in_depth = deck["absorption"] == "in-depth"
    solids = map(comps) do c
        name = Symbol(c["id"])
        ρ = Float64(c["rho"])
        κ = in_depth ? Float64(get(c, "kappa", 0.0)) : 0.0
        SolidComponent(
            name;
            ρ = ρ,
            φ = PHI,
            c = CallableProperty(positive(law(c["cp"]), "$(c["id"]).cp")),
            k = CallableProperty(positive(law(c["k"]), "$(c["id"]).k")),
            ε = Float64(c["emissivity"]),
            α = κ / ρ,
            λ = LAMBDA,
        )
    end
    if substrate == :kaowool
        solids = (
            solids...,
            kaowool(:kaowool, KAOWOOL.κ / KAOWOOL.ρ),
            kaowool(:kaowool_bulk, 0.0),
        )
    end
    gas = GaseousComponent(:gas; M = GAS_M, c = gas_cp, k = 0.02)
    ids = [c["id"] for c in comps]
    index(id) = findfirst(==(id), ids)
    ig = length(solids) + 1
    rho0 = _initial_condensed_density(deck)
    basis_mode = get(deck, "order_basis", "own")
    reactions = map(enumerate(deck["reactions"])) do (j, r)
        ir = index(r["reactant"])
        n = Float64(r["n"])
        A = Float64(r["A"])
        res = get(r, "residue", "")
        y = Float64(get(r, "residue_yield", 0.0))
        basis = nothing
        if n != 1
            if basis_mode == "init"
                basis = 1.0
            elseif basis_mode == "density"
                # r = A_d ρ^n e^(−E/RT) ⇔ A_d ρ₀^(n−1) against the initial density
                basis = 1.0
                A *= Float64(get(deck, "density_basis_rho", rho0))^(n - 1)
            end
        end
        kw = (; A, E = Float64(r["E"]), h = Float64(r["dH"]), n, basis)
        if res == "" || y == 0
            Reaction(Symbol("r$j"), ir => ig; kw...)
        else
            Reaction(Symbol("r$j"), ir => (index(res), ig); yields = (y, 1 - y), kw...)
        end
    end
    return Material(
        name = Symbol(deck["id"]),
        components = (solids..., gas),
        reactions = Tuple(reactions),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
end

_initial_condensed_density(deck) =
    1 / sum(Float64(c["y0"]) / Float64(c["rho"]) for c in deck["components"] if c["y0"] > 0)

# Kaowool PM board (matl-db NIST gasification apparatus README, Table 4):
# ρ 256, c 1070 J/(kg K); k(T) from 260 °C up, extended to 0.0397 W/(m K) at
# 20 °C as Aalto-I's deck does. The cell under the PMMA is opaque
# (κ = 8e6 1/m) and absorbs the radiation the PMMA transmits. The board beneath
# it is transparent: its measured k(T) already holds its own radiation, and an
# opaque cell emitting at its centre temperature would add about 4σT³Δz of
# conductivity on cells as wide as the board's.
const KAOWOOL = (
    ρ = 256.0,
    c = 1070.0,
    ε = 0.9,
    κ = 8e6,
    kT = [293.15, 533.15, 811.15, 1089.15, 1366.15],
    kv = [0.0397, 0.0576, 0.085, 0.125, 0.183],
)
kaowool(name, α) = SolidComponent(
    name;
    ρ = KAOWOOL.ρ,
    φ = 0.8,
    c = KAOWOOL.c,
    k = CallableProperty(T -> _table(KAOWOOL.kT, KAOWOOL.kv, T)),
    ε = KAOWOOL.ε,
    α = α,
    λ = LAMBDA,
)

# Initial concentrations: the virgin components at their mass fractions of
# a mixture of bulk density 1/Σ(Yᵢ/ρᵢ) (volume additive), the rest zero.
# With the substrate the sample lies above `z_back`, the opaque Kaowool cell
# between `z_skin` and `z_back`, and the transparent board below `z_skin`.
function initial_concentrations(deck; z_back = 0.0, z_skin = 0.0, substrate = :none)
    ρmix = _initial_condensed_density(deck)
    ξ = [Float64(c["y0"]) * ρmix for c in deck["components"]]
    fields = Any[(z -> (z > z_back ? x : 0.0)) for x in ξ]
    if substrate == :kaowool
        push!(
            fields,
            z -> (z_skin < z < z_back ? KAOWOOL.ρ : 0.0),
            z -> (z < z_skin ? KAOWOOL.ρ : 0.0),
        )
    end
    return [fields; [z -> 0.0]]
end

radiation_model(deck) =
    deck["absorption"] == "in-depth" ? RADIATIVE_TRANSFER : SURFACE_ABSORPTION

# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------

const SOLVE_KWARGS = (; use_ale = true, deactivation = true)

"""
    gasification_problem(deck, q_kW, L_mm; cells_per_mm = 20, t_end, save_every = 1.0,
                         distribution = Diffuse())

The MaCFP-2 1D gasification target: `q_kW` incident, entering at the surface
mixture emissivity (Kirchhoff): absorbed at the face for a surface deck,
crossing it into the depth for an in-depth deck, diffuse unless
`distribution` is `Collimated()` or `Cone(θ₁, θ₂)`; re-radiation to 293 K; no
convection; perfectly insulated, impermeable back, which reflects the
radiation reaching it; 293 K initially.
"""
function gasification_problem(
    deck,
    q_kW,
    L_mm;
    cells_per_mm = 20,
    t_end,
    save_every = 1.0,
    h = 0.0,
    T_amb = T_AMBIENT,
    distribution = Diffuse(),
)
    L = L_mm * 1e-3
    n = round(Int, cells_per_mm * L_mm)
    nodes = [(i - 1) * L / n for i = 1:(n+1)]
    top = RadiativeFluxBC(q_kW * 1e3; distribution) + RadiativeBC(T_inf = T_amb)
    h > 0 && (top = top + ConvectiveBC(h = h, T_inf = T_amb))
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = AdiabaticBC()),
        mass = (
            top = AmbientGasBC(GasReservoir(P = 0.0, T = T_amb)),
            bottom = ImpermeableBC(),
        ),
    )
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = material(deck),
        bc_set = bcs,
        T_initial = z -> T_amb,
        ξ_initial = initial_concentrations(deck),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:save_every:t_end),
    )
end

"""
    tga_problem(deck, β_Kmin; L = 5e-5, n = 1)

The MaCFP-2 TGA target, 300 → 1000 K at `β_Kmin`: one 50 µm cell whose faces
follow T = 300 + βt; the cell lags its faces by about β(L/2)²/2α, a few mK.
"""
function tga_problem(deck, β_Kmin; L = 5e-5, n = 1, save_dT = 0.5)
    deck = merge(deck, Dict("absorption" => "surface"))
    β = β_Kmin / 60
    t_end = 700 / β
    ramp = DirichletThermalBC(t -> 300.0 + β * t)
    bcs = BoundaryConditionSet(
        thermal = (top = ramp, bottom = ramp),
        mass = (
            top = AmbientGasBC(GasReservoir(P = 0.0, T = T_AMBIENT)),
            bottom = AmbientGasBC(GasReservoir(P = 0.0, T = T_AMBIENT)),
        ),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(L, n),
        material = material(deck),
        bc_set = bcs,
        T_initial = z -> 300.0,
        ξ_initial = initial_concentrations(deck),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:(save_dT/β):t_end),
    )
end

# ---------------------------------------------------------------------------
# MaCFP-3: the NIST gasification apparatus (matl-db Validation_Data README)
# ---------------------------------------------------------------------------

const KAOWOOL_THICKNESS = 0.02852          # 5.72 mm bonded disc + 22.8 mm of four discs
const DISC_AREA = π * (0.0699 / 2)^2       # measured diameter 6.99 cm
# Shared apparatus deck: top h and gas temperature, wall temperature for
# re-radiation (walls 14–18 °C, gas at sample height 20–30 °C, no h given).
const APPARATUS = (h = 10.0, T_gas = 298.15, T_wall = 289.15, T0_mass = 294.05)
# The heater as the sample sees it: uniform intensity from the normal out to 50°.
# A coaxial disc fitted to the README's radial flux profile at 50 kW/m² is seen
# out to 33–62° within the profile's error bars; at 25 kW/m² no disc fits them,
# and 50° serves both fluxes. Diffuse() and Collimated() bound it.
const APPARATUS_HEATER = Cone(0.0, deg2rad(50))

# The experiments (README Table 3): file stem, incident flux, measured PMMA mass [g].
# Mass tests see the flux averaged over the disc (49.3 / 24.7 kW/m²); the back
# thermocouples sit within 0.6 cm of the centre, where it is 50 / 25 kW/m².
const APPARATUS_TESTS = [
    (test = "q50_Mass_R3", q = 49.3, m = 26.689),
    (test = "q50_Mass_R4", q = 49.3, m = 25.643),
    (test = "q50_Mass_R5", q = 49.3, m = 28.113),
    (test = "q50_Temp_R1", q = 50.0, m = 24.303),
    (test = "q50_Temp_R2", q = 50.0, m = 25.136),
    (test = "q50_Temp_R3", q = 50.0, m = 27.877),
    (test = "q50_Temp_R4", q = 50.0, m = 27.499),
    (test = "q50_Temp_R5", q = 50.0, m = 28.187),
    (test = "q25_Temp_R1", q = 25.0, m = 28.093),
    (test = "q25_Temp_R2", q = 25.0, m = 27.750),
    (test = "q25_Temp_R3", q = 25.0, m = 26.077),
]

experiment_path(test) = joinpath(
    matl_db(),
    "PMMA",
    "Validation_Data",
    "NIST_Gasification_Apparatus",
    "MaCFP-PMMA_Gasification_$(test).csv",
)

"""The experiment of `test`: time and the mean of its columns after the first (mass, or the back thermocouples)."""
function read_experiment(test)
    a = read_matl_csv(experiment_path(test))
    y = [
        let v = filter(isfinite, a[i, 2:end])
            isempty(v) ? NaN : sum(v) / length(v)
        end for i in axes(a, 1)
    ]
    return a[:, 1], y
end

"""
    apparatus_problem(deck, test; cells_per_mm = 20, t_end, distribution = APPARATUS_HEATER)

One NIST gasification apparatus test: the PMMA, at the thickness that holds the
test's measured mass at the deck's density on the 6.99 cm disc, bonded to
28.52 mm of Kaowool PM; the test's incident flux from the heater
`distribution` (absorbed at the surface emissivity for a surface deck),
convection (h = 10) to gas at 298 K and re-radiation to walls at 289 K;
adiabatic beneath the Kaowool (its 500 s penetration is about 10 mm). Epoxy
ignored. Initial temperature: the test's first back-thermocouple reading, or
294.05 K (20.9 °C) for the mass tests. Returns the problem and the index of the
node at the PMMA–Kaowool interface.
"""
function apparatus_problem(
    deck,
    t;
    cells_per_mm = 20,
    t_end = t.q < 30 ? 1200.0 : 600.0,
    save_every = 1.0,
    h = APPARATUS.h,
    distribution = APPARATUS_HEATER,
)
    ρmix = _initial_condensed_density(deck)
    L = t.m * 1e-3 / (ρmix * DISC_AREA)
    n = round(Int, cells_per_mm * L * 1e3)
    dz = L / n
    # Kaowool cells grow geometrically (×1.1) downward from the PMMA pitch.
    kw = Float64[0.0]
    w = dz
    while kw[end] < KAOWOOL_THICKNESS
        push!(kw, kw[end] + w)
        w *= 1.1
    end
    kw = kw .* (KAOWOOL_THICKNESS / kw[end])
    depth = reverse(kw)                                   # depths below the interface
    nodes = vcat(KAOWOOL_THICKNESS .- depth, [KAOWOOL_THICKNESS + i * dz for i = 1:n])
    z_skin = nodes[length(kw)-1]                          # lower face of the Kaowool cell under the PMMA
    T0 =
        occursin("Temp", t.test) ? first(filter(isfinite, read_experiment(t.test)[2])) :
        APPARATUS.T0_mass
    top =
        RadiativeFluxBC(t.q * 1e3; distribution) +
        RadiativeBC(T_inf = APPARATUS.T_wall) +
        ConvectiveBC(h = h, T_inf = APPARATUS.T_gas)
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = AdiabaticBC()),
        mass = (
            top = AmbientGasBC(GasReservoir(P = 0.0, T = APPARATUS.T_gas)),
            bottom = ImpermeableBC(),
        ),
    )
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = material(deck; substrate = :kaowool),
        bc_set = bcs,
        T_initial = z -> T0,
        ξ_initial = initial_concentrations(
            deck;
            z_back = KAOWOOL_THICKNESS,
            z_skin,
            substrate = :kaowool,
        ),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:save_every:t_end),
    ),
    length(kw)    # the interface node's index
end

# The inert tests of the same apparatus at 50 kW/m² (README, "Heating of Inert
# Materials"): five Kaowool PM discs painted black (ε 0.95) with thermocouples
# 5.72, 11.44 and 17.16 mm below the top; a 3.175 mm copper disc painted black
# (ε 0.92) on four 5.72 mm discs and one 3.175 mm disc of Kaowool PM, with a
# thermocouple at the copper's back face. Copper: ρ 8933, c 385 J/(kg K),
# k 401 W/(m K).
const INERT_TESTS = [
    (
        test = "Black-Insulation",
        top = nothing,
        ε = 0.95,
        kaowool = 5 * 5.72e-3,
        depths = (5.72e-3, 11.44e-3, 17.16e-3),
        t_end = 1200.0,
    ),
    (
        test = "Black-Copper",
        top = (ρ = 8933.0, c = 385.0, k = 401.0, L = 3.175e-3),
        ε = 0.92,
        kaowool = 4 * 5.72e-3 + 3.175e-3,
        depths = (3.175e-3,),
        t_end = 100.0,
    ),
]

inert_experiment_path(test) = joinpath(
    matl_db(),
    "PMMA",
    "Validation_Data",
    "NIST_Gasification_Apparatus",
    "$(test)_q50_Temp.csv",
)

"""
    inert_problem(t; dz_top = 5e-5)

One inert test on the shared apparatus deck: 50 kW/m² absorbed at the painted
face, convection (h = 10) to gas at 298 K, re-radiation to walls at 289 K,
adiabatic beneath the Kaowool, the experiment's mean first reading as the
initial temperature. The Kaowool cells grow ×1.1 downward from `dz_top`.
Run it under `SURFACE_ABSORPTION`: the painted face is opaque.
"""
function inert_problem(t; dz_top = 5e-5)
    a = read_matl_csv(inert_experiment_path(t.test))
    T0 = sum(a[1, 2:(1+length(t.depths))]) / length(t.depths)
    kw = Float64[0.0]
    w = dz_top
    while kw[end] < t.kaowool
        push!(kw, kw[end] + w)
        w *= 1.1
    end
    kw = kw .* (t.kaowool / kw[end])
    nodes = t.kaowool .- reverse(kw)
    gas = GaseousComponent(:gas; M = GAS_M, c = 0.0, k = 0.02)
    board = SolidComponent(
        :kaowool_bulk;
        ρ = KAOWOOL.ρ,
        φ = 0.8,
        c = KAOWOOL.c,
        k = CallableProperty(T -> _table(KAOWOOL.kT, KAOWOOL.kv, T)),
        ε = t.ε,
        λ = LAMBDA,
    )
    if t.top === nothing
        comps = (board, gas)
        ξ = [z -> KAOWOOL.ρ, z -> 0.0]
    else
        n = 4
        nodes = vcat(nodes, [t.kaowool + i * t.top.L / n for i = 1:n])
        plate = SolidComponent(
            :copper;
            ρ = t.top.ρ,
            φ = PHI,
            c = t.top.c,
            k = t.top.k,
            ε = t.ε,
            λ = LAMBDA,
        )
        comps = (plate, board, gas)
        ξ = [
            z -> (z > t.kaowool ? t.top.ρ : 0.0),
            z -> (z > t.kaowool ? 0.0 : KAOWOOL.ρ),
            z -> 0.0,
        ]
    end
    mat = Material(
        name = Symbol(t.test),
        components = comps,
        reactions = (),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
    top =
        RadiativeFluxBC(50e3) +
        RadiativeBC(T_inf = APPARATUS.T_wall) +
        ConvectiveBC(h = APPARATUS.h, T_inf = APPARATUS.T_gas)
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = AdiabaticBC()),
        mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = mat,
        bc_set = bcs,
        T_initial = z -> T0,
        ξ_initial = ξ,
        tspan = (0.0, t.t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:1.0:t.t_end),
    )
end

"""The temperatures of an inert test at its thermocouple depths below the top."""
function inert_series(sol, t)
    cols = [
        Symbol("T_", replace(@sprintf("%.2f", 1e3d), "." => "_"), "mm_K") for d in t.depths
    ]
    vals = [[temperature_at_depth(sol, k, d) for k in eachindex(sol.t)] for d in t.depths]
    return NamedTuple{(:time_s, :T_top_K, cols...)}((sol.t, sol.surface_T, vals...))
end

"""Apparatus quantities: PMMA mass on the disc [g], MLR, surface T, T at the PMMA–Kaowool interface."""
function apparatus_series(sol, i_if)
    m = condensed_mass(sol)
    mlr = -derivative(sol.t, m) .* 1e3
    mlr[spans_event(sol)] .= NaN
    T_if = map(eachindex(sol.t)) do k
        a, b = sol.T[i_if-1, k], (i_if <= size(sol.T, 1) ? sol.T[i_if, k] : NaN)
        isfinite(b) ? (a + b) / 2 : sol.surface_T[k]
    end
    return (
        time_s = sol.t,
        mass_g = m .* DISC_AREA .* 1e3,
        mlr_g_m2s = mlr,
        T_top_K = sol.surface_T,
        T_back_K = T_if,
    )
end

# ---------------------------------------------------------------------------
# Outputs in the MaCFP quantities
# ---------------------------------------------------------------------------

_condensed(sol) = findall(
    c -> !(c isa GaseousComponent) && !(c.name in (:kaowool, :kaowool_bulk)),
    collect(sol.problem.material.components),
)

function _area(sol)
    a = sol.extras.cross_section_area_history
    return isempty(a) ? ones(length(sol.t)) : a
end

# Condensed mass per unit area [kg/m²] at every saved time.
function condensed_mass(sol)
    ic = _condensed(sol)
    a = _area(sol)
    return [
        sum(sol.m[i, j, k] for i in ic, j in axes(sol.m, 2)) / a[k] for
        k in eachindex(sol.t)
    ]
end

# Temperature `depth` below the moving surface, interpolated over the live
# cell centres (the outermost cell centre for depths shallower than it).
function temperature_at_depth(sol, k, depth)
    zc = view(sol.cell_centers, :, k)
    T = view(sol.T, :, k)
    live = findall(isfinite, zc)
    isempty(live) && return NaN
    ztop = isempty(sol.z) ? maximum(zc[live]) : maximum(filter(isfinite, view(sol.z, :, k)))
    target = ztop - depth
    zl, Tl = zc[live], T[live]
    p = sortperm(zl)
    zl, Tl = zl[p], Tl[p]
    target >= zl[end] && return Tl[end]
    target <= zl[1] && return Tl[1]
    i = searchsortedlast(zl, target)
    w = (target - zl[i]) / (zl[i+1] - zl[i])
    return (1 - w) * Tl[i] + w * Tl[i+1]
end

# Centred differences of a series (one-sided at the ends).
function derivative(t, y)
    d = similar(y)
    n = length(y)
    n < 2 && return fill!(d, NaN)
    d[1] = (y[2] - y[1]) / (t[2] - t[1])
    d[n] = (y[n] - y[n-1]) / (t[n] - t[n-1])
    for i = 2:(n-1)
        d[i] = (y[i+1] - y[i-1]) / (t[i+1] - t[i-1])
    end
    return d
end

# Saved points whose centred difference spans a deactivation event: the event
# vents the layer's last ε_m of mass in one step, which is not a mass loss rate.
function spans_event(sol)
    ev = sol.extras.deactivation_times
    t = sol.t
    lo(i) = t[max(i - 1, 1)]
    hi(i) = t[min(i + 1, length(t))]
    return [any(e -> lo(i) <= e <= hi(i), ev) for i in eachindex(t)]
end

"""The MaCFP-2 gasification quantities of a solve (sample 10 cm × 10 cm)."""
function gasification_series(sol)
    m = condensed_mass(sol)                       # kg/m²
    mlr = -derivative(sol.t, m) .* 1e3            # g/(m² s)
    mlr[spans_event(sol)] .= NaN
    K = eachindex(sol.t)
    t1 = [temperature_at_depth(sol, k, 1e-3) for k in K]
    # UMD's "Top Temperature (1mm)" is the mean over the top 1 mm; BUW-FZJ's
    # front temperature is the mean of devices at 0, 0.5 and 1 mm.
    depth_profile(k, d) = d == 0 ? sol.surface_T[k] : temperature_at_depth(sol, k, d)
    avg1 = [
        sum(depth_profile(k, d) for d in range(0, 1e-3; length = 41)[2:(end-1)]) / 39 *
        (39 / 40) + (depth_profile(k, 0.0) + depth_profile(k, 1e-3)) / 80 for k in K
    ]
    avg3 =
        [(sol.surface_T[k] + temperature_at_depth(sol, k, 0.5e-3) + t1[k]) / 3 for k in K]
    return (
        time_s = sol.t,
        mass_g = m .* SAMPLE_AREA .* 1e3,
        mlr_g_m2s = mlr,
        T_top_K = sol.surface_T,
        T_back_K = sol.extras.back_temperature,
        T_1mm_K = t1,
        thickness_mm = sol.extras.thickness_history .* 1e3,
        T_avg1mm_K = avg1,
        T_3pt_K = avg3,
    )
end

"""The TGA quantities: temperature, mass relative to the initial, its rate per second."""
function tga_series(sol, β_Kmin)
    ev = sol.extras.deactivation_times
    keep = isempty(ev) ? eachindex(sol.t) : findall(<(minimum(ev)), sol.t)
    t = sol.t[keep]
    m = condensed_mass(sol)[keep]
    rel = m ./ m[1]
    return (
        time_s = t,
        T_K = 300.0 .+ (β_Kmin / 60) .* t,
        mass_rel = rel,
        mlr_per_s = -derivative(t, rel),
    )
end

function write_series(path, s::NamedTuple)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(keys(s)), ","))
        for i in eachindex(first(s))
            println(io, join((@sprintf("%.9g", getfield(s, k)[i]) for k in keys(s)), ","))
        end
    end
    return path
end

function read_series(path)
    lines = readlines(path)
    names = Symbol.(split(lines[1], ","))
    rows = [parse.(Float64, split(l, ",")) for l in lines[2:end] if !isempty(strip(l))]
    return NamedTuple{Tuple(names)}(Tuple([r[j] for r in rows] for j in eachindex(names)))
end

# ---------------------------------------------------------------------------
# Comparison: every code's runs against the submissions, the experiments and
# each other
# ---------------------------------------------------------------------------

const RESULTS = joinpath(HERE, "results")

"""A code's saved series of `case`, or `nothing` (`results/<code>/<set>/<variant>/<case>.csv`)."""
function load_run(code, set, variant, case; root = RESULTS)
    p = joinpath(root, code, set, variant, "$case.csv")
    isfile(p) || return nothing
    s = read_series(p)
    length(s.time_s) < 2 && return nothing
    return s
end

"""The submission of a gasification case as the same quantities (its "top" may be 1 mm deep)."""
function load_submission(deck, q, L)
    haskey(deck, "predictions") || return nothing
    p = prediction_path(deck, q, L)
    p === nothing && return nothing
    a = read_matl_csv(p)
    size(a, 1) < 2 && return nothing
    t, m = a[:, 1], a[:, 2]
    mlr = -derivative(t, m) ./ 0.01           # g/s per 0.01 m² → g/(m² s)
    return (time_s = t, mass_g = m, mlr_g_m2s = mlr, T_back_K = a[:, 3], T_top_K = a[:, 4])
end

"""
Metrics of a mass history: peak MLR and its time, the first saved times at or
below 50 % and 10 % of the initial mass, and the final mass fraction.
"""
function metrics(s)
    t, m, mlr = s.time_s, s.mass_g, s.mlr_g_m2s
    ok = findall(isfinite, mlr)
    k = ok[argmax(mlr[ok])]
    m0 = m[1]
    lost(f) = (i = findfirst(x -> x <= (1 - f) * m0, m); i === nothing ? NaN : t[i])
    return (
        peak_mlr = mlr[k],
        t_peak = t[k],
        t50 = lost(0.5),
        t90 = lost(0.9),
        m_end = m[end] / m0,
    )
end

# The submission's "top" temperature: the surface, the mean over the top 1 mm
# (UMD) or the mean of 0, 0.5 and 1 mm (BUW-FZJ).
function top_column(deck)
    k = get(deck, "top_temperature", "surface")
    return k == "avg_1mm" ? :T_avg1mm_K : k == "avg_3pt" ? :T_3pt_K : :T_top_K
end

"""First saved time the mass falls to 0.5 % of its initial value (the series end if never)."""
function burnout_time(s)
    i = findfirst(x -> x <= 0.005 * s.mass_g[1], s.mass_g)
    return i === nothing ? s.time_s[end] : s.time_s[i]
end

"""The time the mass first falls to `f` of its initial value, interpolated between saved times (NaN if never)."""
function crossing(s, f)
    m, t = s.mass_g, s.time_s
    i = findfirst(x -> x <= f * m[1], m)
    (i === nothing || i == 1) && return NaN
    return t[i-1] + (t[i] - t[i-1]) * (m[i-1] - f * m[1]) / (m[i-1] - m[i])
end

function interp(t, y, x)
    (x < t[1] || x > t[end]) && return NaN
    i = clamp(searchsortedlast(t, x), 1, length(t) - 1)
    w = (x - t[i]) / (t[i+1] - t[i])
    return (1 - w) * y[i] + w * y[i+1]
end

# MaCFP-3 quantities: (file key, repeats, column, axis label, title).
const APPARATUS_PANELS = [
    (
        "q50_Mass",
        ["R3", "R4", "R5"],
        :mass_g,
        "mass (g)",
        "49.3 kW/m² (disc average): mass",
    ),
    (
        "q50_MLR",
        ["R3", "R4", "R5"],
        :mlr_g_m2s,
        "MLR (g/m²·s)",
        "49.3 kW/m²: mass loss rate",
    ),
    (
        "q50_Temp",
        ["R1", "R2", "R3", "R4", "R5"],
        :T_back_K,
        "T back (K)",
        "50 kW/m² (centre): back temperature",
    ),
    (
        "q25_Temp",
        ["R1", "R2", "R3"],
        :T_back_K,
        "T back (K)",
        "25 kW/m² (centre): back temperature",
    ),
]

"""One repeat of an apparatus measurement: time and value (mass, MLR or mean back temperature), finite points only."""
function apparatus_experiment(key, r)
    if key == "q50_MLR"
        a = read_matl_csv(
            joinpath(
                matl_db(),
                "PMMA",
                "Validation_Data",
                "NIST_Gasification_Apparatus",
                "MaCFP-PMMA_Gasification_q50_MLR_$(r).csv",
            ),
        )
        t, y = a[:, 1], a[:, 2]
    else
        t, y = read_experiment("$(key)_$r")
    end
    ok = isfinite.(y)
    return t[ok], y[ok]
end

"""The group's own MaCFP-3 prediction files for a quantity ("q50_Mass" etc.)."""
function group_predictions(deck, key)
    v = get(deck, "validation", nothing)
    v === nothing && return String[]
    dir = joinpath(matl_db(), "PMMA", "Validation_Results", dirname(v))
    stem = basename(v) * "_Gasification_" * key
    # The inert calibration tests (Copper, Kaowool) share the q50_Temp stem.
    return sort([
        joinpath(dir, f) for f in readdir(dir) if
        startswith(f, stem) && endswith(f, ".csv") && !occursin(r"Copper|Kaowool", f)
    ])
end

function mean_curve(curves, tgrid)
    return [
        let v = filter(isfinite, [interp(t, y, x) for (t, y) in curves])
            isempty(v) ? NaN : sum(v) / length(v)
        end for x in tgrid
    ]
end

"""
MaCFP's score (matl-db Scripts/NIST_Gasification_Validation.py): the mean
prediction interpolated onto the experiment's times, NRMS = RMS / (max − min)
of the experimental mean, over the shortest repeat.
"""
function nrms(model_curves, exp_curves)
    tgrid = collect(0.0:1.0:minimum(maximum(t) for (t, _) in exp_curves))
    e = mean_curve(exp_curves, tgrid)
    m = mean_curve(model_curves, tgrid)
    ok = isfinite.(e) .& isfinite.(m)
    count(ok) < 10 && return NaN
    return sqrt(sum(abs2, (m .- e)[ok]) / count(ok)) / (maximum(e[ok]) - minimum(e[ok]))
end

"""NRMS of a Pyrolysis.jl apparatus variant (`nothing`: the group's own predictions) for one panel."""
function apparatus_nrms(deck, panel; variant = "base", root = RESULTS)
    key, reps, col = panel
    exp_curves = [apparatus_experiment(key, r) for r in reps]
    if variant === nothing
        curves = [
            (a[:, 1], a[:, 2]) for
            a in (read_matl_csv(f) for f in group_predictions(deck, key)) if size(a, 1) > 1
        ]
    else
        curves = Tuple{Vector{Float64},Vector{Float64}}[]
        for r in reps
            s = load_run(
                "pyrolysis",
                deck["id"],
                variant,
                "apparatus_" * replace(key, "MLR" => "Mass") * "_$r";
                root,
            )
            s === nothing || push!(curves, (s.time_s, getfield(s, col)))
        end
    end
    return isempty(curves) ? NaN : nrms(curves, exp_curves)
end

# ---------------------------------------------------------------------------
# Cost at matched accuracy (a timing campaign of run_timing.py)
# ---------------------------------------------------------------------------

# A run is within tolerance of its code's own finest run when peak MLR and t50
# are within 1 % and the surface and back temperature histories within 2 K; the
# combined error e is the largest of the four ratios (e ≤ 1: within tolerance).
const COST_TOL = (peak = 0.01, t50 = 0.01, T = 2.0)
# Each external code's finest level, its own reference.
const REFLEVEL =
    Dict("fds" => "base_c2.5e-5", "gpyro" => "base_half", "thermakin" => "base_half")

function _read_rows(path)
    lines = filter(!isempty, readlines(path))
    head = split(lines[1], ",")
    return [Dict(zip(head, split(l, ","))) for l in lines[2:end]]
end

"""Combined error of `s` against `ref` over 0 < t ≤ the reference's 90 % mass loss."""
function combined_error(s, ref)
    tb = let b = crossing(ref, 0.1)
        isnan(b) ? ref.time_s[end] : b
    end
    grid = filter(t -> 0 < t <= tb, ref.time_s)
    function maxdiff(col)
        d = filter(
            isfinite,
            [
                abs(
                    interp(s.time_s, getfield(s, col), x) -
                    interp(ref.time_s, getfield(ref, col), x),
                ) for x in grid
            ],
        )
        return isempty(d) ? NaN : maximum(d)
    end
    pk(x) = maximum(filter(isfinite, x.mlr_g_m2s))
    e = (
        abs(pk(s) / pk(ref) - 1) / COST_TOL.peak,
        abs(crossing(s, 0.5) / crossing(ref, 0.5) - 1) / COST_TOL.t50,
        maxdiff(:T_top_K) / COST_TOL.T,
        maxdiff(:T_back_K) / COST_TOL.T,
    )
    return maximum(filter(isfinite, collect(e)); init = 0.0)
end

"""
Every timed run of a campaign in `dir`: (code, set, case, setting, median time
[s], combined error). Pyrolysis.jl's rungs carry the errors time_pyrolysis.jl
scored against its own reference; an external code's levels are scored here
against its finest level, which is timed but never counts as within tolerance.
"""
function cost_rows(dir; root = RESULTS)
    rows = []
    for r in _read_rows(joinpath(dir, "pyrolysis_timing.csv"))
        r["retcode"] == "ERROR" && continue
        e = maximum(
            filter(
                isfinite,
                [
                    parse(Float64, r["err_peak_mlr"]) / COST_TOL.peak,
                    parse(Float64, r["err_t50"]) / COST_TOL.t50,
                    parse(Float64, r["err_T_top_K"]) / COST_TOL.T,
                    parse(Float64, r["err_T_back_K"]) / COST_TOL.T,
                ],
            );
            init = 0.0,
        )
        push!(
            rows,
            (
                "Pyrolysis.jl",
                r["set"],
                r["case"],
                "$(r["cells_per_mm"])/mm ×$(r["tol_factor"])",
                parse(Float64, r["median_s"]),
                e,
            ),
        )
    end
    ext = joinpath(dir, "external_times.csv")
    isfile(ext) || return rows
    for r in _read_rows(ext)
        code, set, level, case = r["code"], r["set"], r["level"], r["case"]
        f(l) = joinpath(root, code, set, l, "gasification_$case.csv")
        isfile(f(level)) && isfile(f(REFLEVEL[code])) || continue
        e =
            level == REFLEVEL[code] ? Inf :
            combined_error(read_series(f(level)), read_series(f(REFLEVEL[code])))
        push!(
            rows,
            (code, set, "gasification_$case", level, parse(Float64, r["median_s"]), e),
        )
    end
    return rows
end

"""
A fixed-step code's time scaled to its own burnout: FDS and ThermaKin step to
the case's end time, Gpyro and Pyrolysis.jl stop at burnout.
"""
function span_scaled(code, set, case, level, t; root = RESULTS)
    code in ("fds", "thermakin") || return t
    tend = Dict("gasification_$(q)kW_$(L)mm" => te for (q, L, te) in GASIFICATION_CASES)
    f = joinpath(root, code, set, level, "$case.csv")
    isfile(f) || return NaN
    b = crossing(read_series(f), 0.005)
    return isnan(b) ? t : t * min(1.0, b / tend[case])
end

# ---------------------------------------------------------------------------
# matl-db prediction files (headers vary by institution: name and unit rows in
# either order, quoted names, BOM, trailing commas). Rows whose first field
# does not parse as a number are headers.
# ---------------------------------------------------------------------------

function read_matl_csv(path)
    rows = Vector{Vector{Float64}}()
    for raw in eachline(path)
        line = replace(raw, '﻿' => "")
        fields = strip.(split(line, r"[,;\t]"))
        isempty(fields) && continue
        v = tryparse(Float64, replace(fields[1], "\"" => ""))
        v === nothing && continue
        push!(
            rows,
            [something(tryparse(Float64, replace(f, "\"" => "")), NaN) for f in fields],
        )
    end
    isempty(rows) && return zeros(0, 0)
    w = minimum(length, rows)
    return reduce(vcat, (permutedims(r[1:w]) for r in rows))
end

function prediction_path(deck, q, L)
    dir = joinpath(matl_db(), "PMMA", "Computational_Results", "Gasification_Predictions")
    p = deck["predictions"]
    for kw in ("kW", "kw")
        f = joinpath(dir, "$(p)_Gasification_$(q)$(kw)_$(L)mm.csv")
        isfile(f) && return f
    end
    return nothing
end

function tga_prediction_path(deck, β)
    preds = get(deck, "tga_predictions", Dict())
    f = get(preds, string(β), nothing)
    f === nothing && return nothing
    p = joinpath(matl_db(), "PMMA", "Computational_Results", "TGA_Predictions", f)
    return isfile(p) ? p : nothing
end

end # module
