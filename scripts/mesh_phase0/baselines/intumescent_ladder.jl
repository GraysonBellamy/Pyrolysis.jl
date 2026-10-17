# Intumescent grid ladder: the performance-audit `rb_intumescent` deck
# (coating → 0.3 foam + 0.7 gas, converted cells swell by 0.3 ρ_coating/ρ_foam)
# at expansion factors 3.6×, 20× and 50× (ρ_foam = 100, 18, 7.2 kg m⁻³) on
# uniform meshes (3.6×: n ∈ {25, 50, 100, 200}, plus 400 when the 200-cell run
# is under three minutes; 20× and 50×: n ∈ {25, 50, 100, 200, 400, 800}) and on
# top-refined stretched meshes (ratio 1.05, 1.10 at n ∈ {50, 100}), 900 s each.
# The 50× variant is also run at n ∈ {100, 200} until the coating is fully
# converted (remaining coating mass below 1e-3 of the initial) or 2700 s,
# whichever comes first.
#
# The stretch-1.10, 100-cell decks throw inside the first step. A probe at the
# end re-runs them with two recording overrides (the residual entry point and
# the cross-section update) and reports, call by call, the time and the
# column-averaged χ̄ every residual evaluation saw until the lateral law threw,
# together with whether any step had been accepted.
#
# Every solve has a wall-clock guard and `maxiters`.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/intumescent_ladder.jl
# With `--analyze-only` the runs are read back from the CSVs and only the
# markdown is regenerated.
#
# Outputs (scripts/mesh_phase0/results/): intumescent_ladder.csv (one row per
# run, written incrementally), intumescent_ladder_failure.csv (one row per
# probed failing deck) and intumescent_ladder.md (tables and findings).

using Pyrolysis
using Logging
using Printf
using Statistics: mean
using OrdinaryDiffEq: DiscreteCallback, SciMLBase
using Pyrolysis.Problem: get_concentration, get_node_positions, chi_range
using Pyrolysis.Residual: cross_section_area_from_state

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESULTS = joinpath(ROOT, "scripts", "mesh_phase0", "results")
mkpath(RESULTS)
const CSV = joinpath(RESULTS, "intumescent_ladder.csv")
const FAIL_CSV = joinpath(RESULTS, "intumescent_ladder_failure.csv")
const MD = joinpath(RESULTS, "intumescent_ladder.md")

const OUTPUT_DT = 4.5
const L0 = 0.003
const T_END = 900.0
const T_LONG = 2700.0
const CONVERTED_FRAC = 1e-3
const RHO_COATING = 1200.0
const FOAM_YIELD = 0.3
const TIMEOUT_S = 600.0
const MAXITERS = 200_000
const N400_BUDGET_S = 180.0
const SOLVE_KW = (
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    min_thickness = 1.5e-4,
    progress = false,
)
const VARIANTS = [
    (label = "3.6x", rho_foam = 100.0),
    (label = "20x", rho_foam = 18.0),
    (label = "50x", rho_foam = 7.2),
]
const UNIFORM_N = Dict(
    "3.6x" => [25, 50, 100, 200],
    "20x" => [25, 50, 100, 200, 400, 800],
    "50x" => [25, 50, 100, 200, 400, 800],
)
const LONG_RUNS = [("50x", 100), ("50x", 200)]
const STRETCHED = [(1.05, 50), (1.05, 100), (1.10, 50), (1.10, 100)]
const FAILURE_PROBE = ("stretched", 100, 1.10)
const T_PROBE = [100.0, 300.0, 600.0, 900.0]

expansion(rho_foam) = FOAM_YIELD * RHO_COATING / rho_foam

# Copy of `rb_intumescent()` (scripts/performance_audit/robustness.jl) with the
# mesh, the foam density and the horizon as arguments; everything else is the
# deck's. Outputs every 4.5 s (201 over 900 s).
function intumescent_problem(mesh; rho_foam, t_end = T_END)
    mat = Material(
        name = :rb_intumescent,
        components = (
            SolidComponent(
                :coating,
                ρ = RHO_COATING,
                c = 1500.0,
                k = 0.3,
                ε = 0.9,
                λ = 1e-6,
            ),
            SolidComponent(:foam, ρ = rho_foam, c = 1200.0, k = 0.1, ε = 0.95, λ = 1e-6),
            GaseousComponent(:gas, M = 0.03, c = 1500.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (
            Reaction(
                :swell,
                1 => (2, 3);
                A = 1e12,
                E = 150e3,
                h = 200e3,
                yields = (FOAM_YIELD, 1 - FOAM_YIELD),
            ),
        ),
        lateral_shrinkage_law = χ -> 1.0 + 0.3 * χ,
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 50e3) +
                      RadiativeBC(T_inf = 300.0) +
                      ConvectiveBC(h = 10.0, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
        ),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:OUTPUT_DT:t_end),
    )
end

make_mesh(kind, n, stretch) =
    kind == "uniform" ? create_uniform_mesh(L0, n, Val(3)) :
    create_stretched_mesh(L0, n, Val(3), stretch; stretch_from = :top)

# Remaining coating mass per unit initial area relative to the initial
# ρ_coating L0, from a state vector (termination test) or from a solution.
function coating_fraction(u, layout, mesh, material)
    z = get_node_positions(u, layout)
    s = 0.0
    @inbounds for i = 1:layout.n_cells
        s += get_concentration(u, layout, i, 1) * (z[i+1] - z[i])
    end
    A = cross_section_area_from_state(mesh, material, u, layout)
    return s * A / (mesh.cross_section_area_initial * RHO_COATING * L0)
end

function coating_fraction_history(sol)
    A = sol.extras.cross_section_area_history
    f = zeros(length(sol.t))
    for k in eachindex(sol.t)
        s = 0.0
        for i = 1:sol.extras.live_cells[k]
            s += sol.ξ[1, i, k] * (sol.z[i+1, k] - sol.z[i, k])
        end
        f[k] = s * A[k] / (A[1] * RHO_COATING * L0)
    end
    return f
end

# Counts warnings (surface-temperature clamp, surface Newton, other), forwards
# warnings and errors, drops the Info noise.
struct WarnCounter{L<:Logging.AbstractLogger} <: Logging.AbstractLogger
    inner::L
    pinned::Base.RefValue{Int}
    newton::Base.RefValue{Int}
    other::Base.RefValue{Int}
end
WarnCounter(inner) = WarnCounter(inner, Ref(0), Ref(0), Ref(0))
Logging.min_enabled_level(::WarnCounter) = Logging.Warn
Logging.shouldlog(::WarnCounter, level, args...) = level >= Logging.Warn
Logging.catch_exceptions(::WarnCounter) = false
function Logging.handle_message(l::WarnCounter, level, message, args...; kwargs...)
    if level == Logging.Warn
        msg = string(message)
        if occursin("pinned at the physical clamp", msg)
            l.pinned[] += 1
        elseif occursin("Surface-temperature Newton did not converge", msg)
            l.newton[] += 1
        else
            l.other[] += 1
        end
        return nothing
    end
    return Logging.handle_message(l.inner, level, message, args...; kwargs...)
end

# Wall-clock guard and optional state-based stop. The condition runs after
# every accepted step, records its time and counts it; it is true only once
# the budget is spent or `stop(u)` holds, so the callback does not touch the
# integrator before that.
function guarded_solve(problem; limit_s = TIMEOUT_S, maxiters = MAXITERS, stop = nothing)
    t0 = time()
    hit = Ref(false)
    converted = Ref(false)
    t_last = Ref(NaN)
    n_accepted = Ref(0)
    guard = DiscreteCallback(
        (u, t, integ) -> begin
            t_last[] = t
            n_accepted[] += 1
            time() - t0 > limit_s || (stop !== nothing && stop(u))
        end,
        integ -> begin
            if time() - t0 > limit_s
                hit[] = true
            else
                converted[] = true
            end
            SciMLBase.terminate!(integ)
        end;
        save_positions = (false, false),
    )
    wc = WarnCounter(Logging.current_logger())
    sol = try
        Logging.with_logger(wc) do
            solve(problem; callback = guard, maxiters, SOLVE_KW...)
        end
    catch err
        (; error = err)
    end
    return (;
        sol,
        wall = time() - t0,
        hit = hit[],
        converted = converted[],
        t_last = t_last[],
        n_accepted = n_accepted[],
        wc,
    )
end

function interp_at(t, y, tq)
    (tq < t[1] || tq > t[end]) && return NaN
    j = searchsortedlast(t, tq)
    j == length(t) && return y[end]
    w = (tq - t[j]) / (t[j+1] - t[j])
    return (1 - w) * y[j] + w * y[j+1]
end

const COLUMNS = [
    "variant",
    "expansion",
    "rho_foam",
    "mesh",
    "stretch",
    "n",
    "horizon",
    "retcode",
    "timed_out",
    "converted",
    "t_reached",
    "t_last_accepted",
    "accepted_seen",
    "wall_s",
    "steps",
    "rejects",
    "nf",
    "njacs",
    "nonliniter",
    "nonlinfail",
    "warn_pinned",
    "warn_newton",
    "warn_other",
    "Ts_100",
    "Ts_300",
    "Ts_600",
    "Ts_900",
    "Tb_900",
    "L_900",
    "L_end",
    "coating_frac_900",
    "coating_frac_end",
    "t_converted",
    "peak_mlr",
    "t_peak",
    "dz0_min",
    "dz0_max",
    "dz_min",
    "dz_max",
    "dz_ratio",
]
const STRING_COLS = ("variant", "mesh", "retcode")
const BOOL_COLS = ("timed_out", "converted")
const INT_COLS = (
    "n",
    "accepted_seen",
    "steps",
    "rejects",
    "nf",
    "njacs",
    "nonliniter",
    "nonlinfail",
    "warn_pinned",
    "warn_newton",
    "warn_other",
)
const QOIS = [:Ts_100, :Ts_300, :Ts_600, :Ts_900, :Tb_900, :L_end, :peak_mlr, :t_peak]
# Differences below these are within the time-integration tolerance (reltol
# 1e-6 on ~900 K states) or one output interval (t_peak): no order is formed.
const FLOOR = Dict(
    :Ts_100 => 1e-2,
    :Ts_300 => 1e-2,
    :Ts_600 => 1e-2,
    :Ts_900 => 1e-2,
    :Tb_900 => 1e-2,
    :L_end => 1e-6,
    :peak_mlr => 1e-2,
    :t_peak => 4.5,
)
const UNIT = Dict(
    :Ts_100 => "K",
    :Ts_300 => "K",
    :Ts_600 => "K",
    :Ts_900 => "K",
    :Tb_900 => "K",
    :L_end => "m",
    :peak_mlr => "g/m²/s",
    :t_peak => "s",
)

function run_one(variant, kind, n, stretch; t_end = T_END, stop_on_conversion = false)
    row = Dict{String,Any}(c => NaN for c in COLUMNS)
    row["variant"] = variant.label
    row["expansion"] = expansion(variant.rho_foam)
    row["rho_foam"] = variant.rho_foam
    row["mesh"] = kind
    row["stretch"] = kind == "uniform" ? 1.0 : stretch
    row["n"] = n
    row["horizon"] = t_end
    row["timed_out"] = false
    row["converted"] = false
    row["accepted_seen"] = 0
    row["retcode"] = "not run"
    try
        mesh = make_mesh(kind, n, stretch)
        dz0 = diff(mesh.face_positions)
        row["dz0_min"], row["dz0_max"] = extrema(dz0)
        problem = intumescent_problem(mesh; rho_foam = variant.rho_foam, t_end)
        stop = nothing
        if stop_on_conversion
            # Layout, initial dry masses and material of an identical problem;
            # the stop test reads only the integrator's state vector.
            p2 = intumescent_problem(
                make_mesh(kind, n, stretch);
                rho_foam = variant.rho_foam,
                t_end,
            )
            def, ws, _ =
                initialize_state(p2; use_ale = true, radiation_model = SURFACE_ABSORPTION)
            layout, mesh2, material = ws.layout, ws.mesh, def.material
            stop = u -> coating_fraction(u, layout, mesh2, material) < CONVERTED_FRAC
        end
        r = guarded_solve(problem; stop)
        row["wall_s"] = r.wall
        row["timed_out"] = r.hit
        row["converted"] = r.converted
        row["t_last_accepted"] = r.t_last
        row["accepted_seen"] = r.n_accepted
        row["warn_pinned"] = r.wc.pinned[]
        row["warn_newton"] = r.wc.newton[]
        row["warn_other"] = r.wc.other[]
        sol = r.sol
        if sol isa NamedTuple
            row["retcode"] = "Exception: " * first(sprint(showerror, sol.error), 160)
            return row
        end
        st = sol.raw_solution.stats
        row["retcode"] = string(sol.retcode)
        row["t_reached"] = sol.t[end]
        row["steps"] = st.naccept
        row["rejects"] = st.nreject
        row["nf"] = st.nf
        row["njacs"] = st.njacs
        row["nonliniter"] = st.nnonliniter
        row["nonlinfail"] = st.nnonlinconvfail
        t = sol.t
        for tq in T_PROBE
            row["Ts_$(Int(tq))"] = interp_at(t, sol.surface_T, tq)
        end
        row["Tb_900"] = interp_at(t, sol.extras.back_temperature, 900.0)
        row["L_900"] = interp_at(t, sol.extras.thickness_history, 900.0)
        row["L_end"] = sol.extras.thickness_history[end]
        f = coating_fraction_history(sol)
        row["coating_frac_900"] = interp_at(t, f, 900.0)
        row["coating_frac_end"] = f[end]
        row["t_converted"] = r.converted ? r.t_last : NaN
        mlr = sol.mass_loss_rate
        i = argmax(mlr)
        row["peak_mlr"] = 1e3 * mlr[i]
        row["t_peak"] = t[i]
        n_live = sol.extras.live_cells[end]
        dz = diff(sol.z[1:(n_live+1), end])
        row["dz_min"], row["dz_max"] = extrema(dz)
        row["dz_ratio"] = row["dz_max"] / row["dz_min"]
    catch err
        row["retcode"] = "Exception (setup): " * first(sprint(showerror, err), 160)
    end
    return row
end

csvcell(x) =
    x isa AbstractString ? "\"" * replace(x, "\"" => "'", "\n" => " ") * "\"" : string(x)
function write_row(io, row, columns = COLUMNS)
    println(io, join((csvcell(row[c]) for c in columns), ","))
    flush(io)
end
function report(row)
    println(
        @sprintf(
            "%-5s %-9s stretch=%.2f n=%4d horizon=%4.0f: %s%s%s t=%s wall=%.1f s steps=%s rej=%s nlfail=%s warn=%s/%s/%s Ts900=%s L_end=%s coating_end=%s dz_ratio=%s",
            row["variant"],
            row["mesh"],
            row["stretch"],
            row["n"],
            row["horizon"],
            first(row["retcode"], 60),
            row["timed_out"] ? " (timed out)" : "",
            row["converted"] ? " (converted)" : "",
            row["t_reached"],
            row["wall_s"],
            row["steps"],
            row["rejects"],
            row["nonlinfail"],
            row["warn_pinned"],
            row["warn_newton"],
            row["warn_other"],
            row["Ts_900"],
            row["L_end"],
            row["coating_frac_end"],
            row["dz_ratio"]
        ),
    )
    flush(stdout)
end

function parse_csv_line(line)
    fields = String[]
    buf = IOBuffer()
    inq = false
    for ch in line
        if ch == '"'
            inq = !inq
        elseif ch == ',' && !inq
            push!(fields, String(take!(buf)))
        else
            write(buf, ch)
        end
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_runs(path)
    lines = readlines(path)
    parse_csv_line(lines[1]) == COLUMNS || error("unexpected header in $path")
    out = Dict{String,Any}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        r = Dict{String,Any}()
        for (c, x) in zip(COLUMNS, parse_csv_line(line))
            r[c] =
                c in STRING_COLS ? x :
                c in BOOL_COLS ? parse(Bool, x) :
                c in INT_COLS ? (x == "NaN" ? NaN : parse(Int, x)) : parse(Float64, x)
        end
        push!(out, r)
    end
    return out
end

# ----------------------------------------------------------------------------
# First-step failure probe: process-local recording overrides
#
# `residual!` records the time and the scalar type (Float64 = residual or
# Newton evaluation, Dual = Jacobian pass) of every call; the WithChi
# cross-section update records the χ̄ it was handed, tagged with the index of
# the residual call it belongs to. Both are copies of the package methods
# (src/Residual/residual.jl, src/Problem/state_cache.jl) with the recording
# lines as the only addition; the error message is the package's.
# ----------------------------------------------------------------------------
const PROBE_T = Float64[]
const PROBE_DUAL = Bool[]
const PROBE_CHI = Float64[]
const PROBE_CHI_CALL = Int[]

function install_failure_probe()
    @eval Pyrolysis.Residual function residual!(
        du::AbstractVector,
        u::AbstractVector,
        def::ProblemDef,
        t::Real,
        ws::Workspace,
    )
        push!(Main.PROBE_T, Float64(t))
        push!(Main.PROBE_DUAL, eltype(u) <: Main.Pyrolysis.Problem.ForwardDiff.Dual)
        _residual!(du, u, def, t, ws, mesh_mode(def.mode))
        return nothing
    end
    @eval Pyrolysis.Problem function update_cross_section_area!(
        cache::StateCache{NC,Tu},
        material,
        mesh,
        u::AbstractVector,
        layout::StateLayout{NC},
        ::WithChi,
    ) where {NC,Tu}
        if material.lateral_shrinkage_law === nothing
            return (cache.A_section[], zero(Tu), 0.0)
        end
        χ̄, M = _chi_bar_from_state(u, layout, mesh)
        A = mesh.cross_section_area_initial * material.lateral_shrinkage_law(χ̄)
        push!(Main.PROBE_CHI, ForwardDiff.value(χ̄))
        push!(Main.PROBE_CHI_CALL, length(Main.PROBE_T))
        # string(): a `$` inside this @eval body would be interpolated at eval time.
        ForwardDiff.value(A) > 0.0 || error(
            string(
                "lateral_shrinkage_law(χ̄ = ",
                ForwardDiff.value(χ̄),
                ") gives a non-positive cross-section area A = ",
                ForwardDiff.value(A),
                " m². The law must stay strictly positive on χ ∈ [0, 1]; ",
                "non-positive A makes every cell volume non-positive.",
            ),
        )
        cache.A_section[] = A
        return (A, χ̄, M)
    end
    return nothing
end

const FAIL_COLUMNS = [
    "variant",
    "mesh",
    "stretch",
    "n",
    "init_area",
    "init_chi_min",
    "init_chi_max",
    "init_residual_finite",
    "retcode",
    "accepted_seen",
    "t_last_accepted",
    "wall_s",
    "residual_calls",
    "jacobian_assemblies",
    "dual_passes",
    "calls_at_t0",
    "t_throw",
    "t_throw_over_dt0",
    "throw_pass",
    "chi_bar_throw",
    "chi_bar_message",
    "first_nonzero_call",
    "t_first_nonzero",
    "chi_bar_first_nonzero",
    "chi_bar_max_abs_before_throw",
    "warn_pinned",
    "warn_newton",
    "warn_other",
    "chi_trace",
]

function run_failure_probe(variant)
    kind, n, stretch = FAILURE_PROBE
    row = Dict{String,Any}(c => NaN for c in FAIL_COLUMNS)
    row["variant"] = variant.label
    row["mesh"] = kind
    row["stretch"] = stretch
    row["n"] = n
    row["chi_trace"] = ""
    # Initialization: area, χ and residual of the initial state (the
    # reviewer's check).
    prob0 = intumescent_problem(make_mesh(kind, n, stretch); rho_foam = variant.rho_foam)
    def, ws, u0 =
        initialize_state(prob0; use_ale = true, radiation_model = SURFACE_ABSORPTION)
    du = similar(u0)
    Pyrolysis.Residual.residual!(du, u0, def, 0.0, ws)
    row["init_area"] = ws.mesh.cross_section_area
    row["init_chi_min"], row["init_chi_max"] = extrema(u0[chi_range(ws.layout)])
    row["init_residual_finite"] = all(isfinite, du)
    # The solve, recorded from its first residual call.
    empty!(PROBE_T)
    empty!(PROBE_DUAL)
    empty!(PROBE_CHI)
    empty!(PROBE_CHI_CALL)
    problem = intumescent_problem(make_mesh(kind, n, stretch); rho_foam = variant.rho_foam)
    r = guarded_solve(problem)
    row["wall_s"] = r.wall
    row["accepted_seen"] = r.n_accepted
    row["t_last_accepted"] = r.t_last
    row["warn_pinned"] = r.wc.pinned[]
    row["warn_newton"] = r.wc.newton[]
    row["warn_other"] = r.wc.other[]
    msg = r.sol isa NamedTuple ? sprint(showerror, r.sol.error) : ""
    row["retcode"] =
        r.sol isa NamedTuple ? "Exception: " * first(msg, 160) : string(r.sol.retcode)
    m = match(r"χ̄ = ([-+0-9.eE]+)", msg)
    row["chi_bar_message"] = m === nothing ? NaN : parse(Float64, m.captures[1])
    row["residual_calls"] = length(PROBE_T)
    # Every ALE residual call updates the area exactly once; a second update
    # tagged with the same call index comes from outside `residual!`: the
    # structured Jacobian assembly (src/Jacobian/structured_backend.jl) calls
    # `update_cross_section_area!` directly at the step's base state.
    dup = [PROBE_CHI_CALL[k] in view(PROBE_CHI_CALL, 1:(k-1)) for k in eachindex(PROBE_CHI)]
    pass(k) =
        dup[k] ? "A" : (PROBE_CHI_CALL[k] >= 1 && PROBE_DUAL[PROBE_CHI_CALL[k]]) ? "J" : "r"
    row["jacobian_assemblies"] = count(dup)
    row["dual_passes"] = count(PROBE_DUAL)
    row["calls_at_t0"] = count(==(0.0), PROBE_T)
    row["throw_pass"] = ""
    if !isempty(PROBE_CHI)
        c_last = PROBE_CHI_CALL[end]
        row["t_throw"] = c_last >= 1 ? PROBE_T[c_last] : NaN
        row["t_throw_over_dt0"] = row["t_throw"] / 1e-3
        row["throw_pass"] = pass(length(PROBE_CHI))
        row["chi_bar_throw"] = PROBE_CHI[end]
        k = findfirst(!=(0.0), PROBE_CHI)
        if k !== nothing
            c = PROBE_CHI_CALL[k]
            row["first_nonzero_call"] = c
            row["t_first_nonzero"] = c >= 1 ? PROBE_T[c] : NaN
            row["chi_bar_first_nonzero"] = PROBE_CHI[k]
        else
            row["first_nonzero_call"] = 0
        end
        row["chi_bar_max_abs_before_throw"] =
            length(PROBE_CHI) > 1 ? maximum(abs, PROBE_CHI[1:(end-1)]) : 0.0
        row["chi_trace"] = join(
            (
                @sprintf(
                    "c%d t=%.4g %s χ̄=%.4g",
                    PROBE_CHI_CALL[k],
                    PROBE_CHI_CALL[k] >= 1 ? PROBE_T[PROBE_CHI_CALL[k]] : NaN,
                    pass(k),
                    PROBE_CHI[k]
                ) for k in eachindex(PROBE_CHI)
            ),
            ";",
        )
    end
    return row
end

function read_failures(path)
    isfile(path) || return Dict{String,Any}[]
    lines = readlines(path)
    parse_csv_line(lines[1]) == FAIL_COLUMNS || error("unexpected header in $path")
    out = Dict{String,Any}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        r = Dict{String,Any}()
        for (c, x) in zip(FAIL_COLUMNS, parse_csv_line(line))
            r[c] =
                c in ("variant", "mesh", "retcode", "chi_trace", "throw_pass") ? x :
                c == "init_residual_finite" ? parse(Bool, x) :
                c in (
                    "n",
                    "accepted_seen",
                    "residual_calls",
                    "jacobian_assemblies",
                    "dual_passes",
                    "calls_at_t0",
                    "first_nonzero_call",
                    "warn_pinned",
                    "warn_newton",
                    "warn_other",
                ) ? (x == "NaN" ? NaN : parse(Int, x)) : parse(Float64, x)
        end
        push!(out, r)
    end
    return out
end

# ----------------------------------------------------------------------------
# Runs
# ----------------------------------------------------------------------------
t_script = time()
rows = Dict{String,Any}[]
frows = Dict{String,Any}[]
if "--analyze-only" in ARGS
    rows = read_runs(CSV)
    frows = read_failures(FAIL_CSV)
    println("read ", length(rows), " runs from ", CSV, " and ", length(frows), " probes")
else
    println("intumescent ladder: timeout $TIMEOUT_S s, maxiters $MAXITERS")
    flush(stdout)
    let r = run_one(VARIANTS[1], "uniform", 25, 1.0)
        println("warmup: ", r["retcode"], " ", @sprintf("%.1f", r["wall_s"]), " s")
        flush(stdout)
    end
    open(CSV, "w") do io
        println(io, join(COLUMNS, ","))
        flush(io)
        for v in VARIANTS
            ns = copy(UNIFORM_N[v.label])
            while !isempty(ns)
                n = popfirst!(ns)
                r = run_one(v, "uniform", n, 1.0)
                push!(rows, r)
                write_row(io, r)
                report(r)
                if v.label == "3.6x" &&
                   n == 200 &&
                   r["retcode"] == "Success" &&
                   r["wall_s"] < N400_BUDGET_S
                    push!(ns, 400)
                end
            end
            for (stretch, n) in STRETCHED
                r = run_one(v, "stretched", n, stretch)
                push!(rows, r)
                write_row(io, r)
                report(r)
            end
            for (label, n) in LONG_RUNS
                label == v.label || continue
                r = run_one(v, "uniform", n, 1.0; t_end = T_LONG, stop_on_conversion = true)
                push!(rows, r)
                write_row(io, r)
                report(r)
            end
        end
    end
end

# ----------------------------------------------------------------------------
# Failure probe (after the ladder: the overrides record on every residual call)
# ----------------------------------------------------------------------------
if !("--analyze-only" in ARGS)
    install_failure_probe()
end
if !("--analyze-only" in ARGS)
    open(FAIL_CSV, "w") do io
        println(io, join(FAIL_COLUMNS, ","))
        flush(io)
        for v in VARIANTS
            r = run_failure_probe(v)
            push!(frows, r)
            write_row(io, r, FAIL_COLUMNS)
            println(
                @sprintf(
                    "probe %-5s stretch=%.2f n=%d: init area=%s χ∈[%s,%s] residual finite=%s | %s | accepted=%d residual calls=%d (Jacobian assemblies %d, dual passes %d, t=0: %d) throw at t=%.4g = %.3g dt0 (pass %s) χ̄=%.4g (message %.4g) first χ̄≠0 at call %s (t=%.4g, χ̄=%.4g) max|χ̄| before throw=%.4g warn=%s/%s/%s wall=%.2f s",
                    r["variant"],
                    r["stretch"],
                    r["n"],
                    r["init_area"],
                    r["init_chi_min"],
                    r["init_chi_max"],
                    r["init_residual_finite"],
                    first(r["retcode"], 70),
                    r["accepted_seen"],
                    r["residual_calls"],
                    r["jacobian_assemblies"],
                    r["dual_passes"],
                    r["calls_at_t0"],
                    r["t_throw"],
                    r["t_throw_over_dt0"],
                    r["throw_pass"],
                    r["chi_bar_throw"],
                    r["chi_bar_message"],
                    r["first_nonzero_call"],
                    r["t_first_nonzero"],
                    r["chi_bar_first_nonzero"],
                    r["chi_bar_max_abs_before_throw"],
                    r["warn_pinned"],
                    r["warn_newton"],
                    r["warn_other"],
                    r["wall_s"]
                ),
            )
            println("  trace: ", r["chi_trace"])
            flush(stdout)
        end
    end
end

# ----------------------------------------------------------------------------
# Convergence across the uniform ladder and grading comparison
# ----------------------------------------------------------------------------
fmtv(x) =
    x isa AbstractString ? x :
    !isfinite(x) ? string(x) :
    (abs(x) >= 1e5 || (0 < abs(x) < 1e-3) ? @sprintf("%.4e", x) : @sprintf("%.6g", x))
ok(r) = r["retcode"] == "Success" && !r["timed_out"]
base_horizon(r) = r["horizon"] == T_END
uniform_rows(v) = sort(
    filter(r -> r["variant"] == v.label && r["mesh"] == "uniform" && base_horizon(r), rows),
    by = r -> r["n"],
)
stretched_rows(v) = sort(
    filter(
        r -> r["variant"] == v.label && r["mesh"] == "stretched" && base_horizon(r),
        rows,
    ),
    by = r -> (r["stretch"], r["n"]),
)
long_rows(v) =
    sort(filter(r -> r["variant"] == v.label && !base_horizon(r), rows), by = r -> r["n"])

# Observed orders from successive differences; NaN when either difference is
# below the QoI's floor.
function observed_orders(vals, floor)
    d = diff(vals)
    return [
        (abs(d[k]) > floor && abs(d[k+1]) > floor) ? log2(abs(d[k]) / abs(d[k+1])) : NaN for
        k = 1:(length(d)-1)
    ]
end

# First n after which every successive difference stays below the floor.
function converged_at(ns, vals, floor)
    d = abs.(diff(vals))
    for k in eachindex(d)
        all(<(floor), d[k:end]) && return ns[k]
    end
    return nothing
end

fmt_or_dash(x) = isnan(x) ? "–" : fmtv(x)

open(MD, "w") do io
    println(io, "# Intumescent grid ladder\n")
    println(
        io,
        "Generated by `scripts/mesh_phase0/baselines/intumescent_ladder.jl`; rows in `intumescent_ladder.csv`, ",
        "failure-probe rows in `intumescent_ladder_failure.csv`.\n",
    )
    println(
        io,
        "Deck: performance-audit `rb_intumescent` (3 mm coating, 50 kW/m² + reradiation + convection, ",
        "adiabatic back, coating → 0.3 foam + 0.7 gas, lateral law 1 + 0.3 χ, outputs every 4.5 s) ",
        "with ρ_foam = 100, 18, 7.2 kg m⁻³ giving volumetric expansion 0.3 ρ_coating/ρ_foam = ",
        join((@sprintf("%.1f", expansion(v.rho_foam)) for v in VARIANTS), ", "),
        ". Uniform ladders: 3.6× at n = ",
        join(UNIFORM_N["3.6x"], "/"),
        " (400 added when the 200-cell run is under $(Int(N400_BUDGET_S)) s), 20× and 50× at n = ",
        join(UNIFORM_N["20x"], "/"),
        "; stretched meshes (finest cells at the exposed surface) at ratio 1.05 and 1.10, n = 50 and 100; ",
        "all over 900 s. Long-horizon runs: 50× at n = 100 and 200 until the remaining coating mass ",
        "(Σ ξ_coating Δz · A/A₀, relative to ρ_coating L₀) falls below $(CONVERTED_FRAC) or t = $(Int(T_LONG)) s. ",
        "Solver: `use_ale = true`, `SURFACE_ABSORPTION`, `min_thickness = 1.5e-4`, defaults otherwise ",
        "(ApproxSolve + KLU, KenCarp4, abstol 1e-8, reltol 1e-6), maxiters $MAXITERS, wall-clock guard $TIMEOUT_S s. ",
        "Surface and back temperatures and the thickness at 900 s are linearly interpolated between outputs; `t_peak` and the peak MLR ",
        "are read on the output grid. `warn` counts surface-temperature clamp / surface-Newton / other warnings. ",
        "Wall times are indicative (shared machine, another Julia process running concurrently).\n",
    )
    println(io, "## Runs\n")
    println(
        io,
        "| variant | mesh | stretch | n | horizon s | retcode | timed out | converted | t reached | wall s | steps | rejects | nlfail | warn | Ts(100) | Ts(300) | Ts(600) | Ts(900) | Tb(900) | L(900) mm | L_end mm | coating left (900 / end) | peak MLR g/m²/s | t_peak | dz0 min/max mm | dz_end min/max mm | dz ratio |",
    )
    println(
        io,
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    )
    for r in rows
        println(
            io,
            @sprintf(
                "| %s | %s | %.2f | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s/%s/%s | %s | %s | %s | %s | %s | %s | %s | %s / %s | %s | %s | %s/%s | %s/%s | %s |",
                r["variant"],
                r["mesh"],
                r["stretch"],
                r["n"],
                r["horizon"],
                first(r["retcode"], 70),
                r["timed_out"],
                r["converted"],
                fmtv(r["t_reached"]),
                fmtv(r["wall_s"]),
                r["steps"],
                r["rejects"],
                r["nonlinfail"],
                r["warn_pinned"],
                r["warn_newton"],
                r["warn_other"],
                fmtv(r["Ts_100"]),
                fmtv(r["Ts_300"]),
                fmtv(r["Ts_600"]),
                fmtv(r["Ts_900"]),
                fmtv(r["Tb_900"]),
                fmtv(1e3 * r["L_900"]),
                fmtv(1e3 * r["L_end"]),
                fmtv(r["coating_frac_900"]),
                fmtv(r["coating_frac_end"]),
                fmtv(r["peak_mlr"]),
                fmtv(r["t_peak"]),
                fmtv(1e3 * r["dz0_min"]),
                fmtv(1e3 * r["dz0_max"]),
                fmtv(1e3 * r["dz_min"]),
                fmtv(1e3 * r["dz_max"]),
                fmtv(r["dz_ratio"])
            )
        )
    end

    # ---- convergence tables per variant ------------------------------------
    println(io, "\n## Convergence across the uniform ladder (900 s)\n")
    println(
        io,
        "Successive differences Δ(n→2n) = q(2n) − q(n) and observed order p = log2(|Δ(n→2n)| / |Δ(2n→4n)|), ",
        "formed only when both differences exceed the QoI's floor. Floors: ",
        join(("$q $(FLOOR[q]) $(UNIT[q])" for q in QOIS), ", "),
        " (temperatures: ten times reltol × 900 K; thickness: 1 µm; t_peak: one output interval). ",
        "`at floor from` is the first n after which every successive difference stays below the floor; ",
        "`moving` means the finest pair still differs by more than the floor.\n",
    )
    conv = Dict{String,Dict{Symbol,Any}}()
    for v in VARIANTS
        ur = filter(ok, uniform_rows(v))
        ns = [r["n"] for r in ur]
        println(
            io,
            "### $(v.label) (ρ_foam = $(v.rho_foam), expansion $(@sprintf("%.1f", expansion(v.rho_foam))))\n",
        )
        conv[v.label] = Dict{Symbol,Any}()
        if length(ns) < 2
            println(io, "Fewer than two successful uniform runs.\n")
            continue
        end
        println(
            io,
            "| QoI | floor | ",
            join(("n = $n" for n in ns), " | "),
            " | ",
            join(("Δ($(ns[k])→$(ns[k+1]))" for k = 1:(length(ns)-1)), " | "),
            " | ",
            join(("p($(ns[k]),$(ns[k+1]),$(ns[k+2]))" for k = 1:(length(ns)-2)), " | "),
            " | status |",
        )
        println(
            io,
            "|---|---|",
            "---|"^(length(ns) + (length(ns) - 1) + (length(ns) - 2) + 1),
        )
        for q in QOIS
            vals = [Float64(r[string(q)]) for r in ur]
            d = diff(vals)
            p = observed_orders(vals, FLOOR[q])
            c = converged_at(ns, vals, FLOOR[q])
            conv[v.label][q] = (; ns, vals, d, p, c)
            println(
                io,
                "| $q | $(FLOOR[q]) $(UNIT[q]) | ",
                join(fmtv.(vals), " | "),
                " | ",
                join(fmtv.(d), " | "),
                " | ",
                join(fmt_or_dash.(p), " | "),
                " | ",
                c === nothing ? "moving" : "at floor from n = $c",
                " |",
            )
        end
        println(io)
    end

    # ---- status at the finest pair, per quantity and expansion factor ------
    println(io, "### Finest pair of each ladder, per quantity\n")
    println(
        io,
        "| QoI | floor | ",
        join(
            (
                "$(v.label): Δ($(conv[v.label][:Ts_100].ns[end-1])→$(conv[v.label][:Ts_100].ns[end])) / p / status"
                for v in VARIANTS if haskey(conv[v.label], :Ts_100)
            ),
            " | ",
        ),
        " |",
    )
    println(io, "|---|---|", "---|"^count(v -> haskey(conv[v.label], :Ts_100), VARIANTS))
    for q in QOIS
        cells = String[]
        for v in VARIANTS
            haskey(conv[v.label], q) || continue
            c = conv[v.label][q]
            push!(
                cells,
                @sprintf(
                    "%s %s / %s / %s",
                    fmtv(c.d[end]),
                    UNIT[q],
                    isempty(c.p) ? "–" : fmt_or_dash(c.p[end]),
                    abs(c.d[end]) < FLOOR[q] ?
                    (c.c === nothing ? "below floor" : "at floor from n = $(c.c)") :
                    "moving"
                )
            )
        end
        println(io, "| $q | $(FLOOR[q]) $(UNIT[q]) | ", join(cells, " | "), " |")
    end

    # ---- long-horizon runs --------------------------------------------------
    println(io, "\n## 50× to full conversion\n")
    lr = long_rows(VARIANTS[3])
    println(
        io,
        "Stop when the remaining coating mass falls below $(CONVERTED_FRAC) of the initial, or at $(Int(T_LONG)) s. ",
        "`t_converted` is the integrator time of the first accepted step at which the test held (NaN when the horizon ended first); ",
        "it is resolved only to the step size there (the runs average about $(Int(round(1926 / 330))) s per step). `t reached` is the last saved ",
        "output time (the solution's output grid stops at the last output before the termination), and `L_end`, `coating left at end` ",
        "are read there, so both runs are compared at the same time. The first long run's wall time includes the compilation of the stop test.\n",
    )
    println(
        io,
        "| n | retcode | converted | t_converted s | t reached | steps | wall s | coating left at 900 s | coating left at end | L(900) mm | L_end mm | Ts(900) | Tb(900) |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in lr
        println(
            io,
            @sprintf(
                "| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                r["n"],
                first(r["retcode"], 70),
                r["converted"],
                fmtv(r["t_converted"]),
                fmtv(r["t_reached"]),
                r["steps"],
                fmtv(r["wall_s"]),
                fmtv(r["coating_frac_900"]),
                fmtv(r["coating_frac_end"]),
                fmtv(1e3 * r["L_900"]),
                fmtv(1e3 * r["L_end"]),
                fmtv(r["Ts_900"]),
                fmtv(r["Tb_900"])
            )
        )
    end
    long_ok = filter(
        ok_r -> ok_r["retcode"] in ("Success", "Terminated") && !ok_r["timed_out"],
        lr,
    )
    if length(long_ok) == 2
        a, b = long_ok
        println(
            io,
            "\nDifference n = $(a["n"]) → $(b["n"]): L_end ",
            fmtv(b["L_end"] - a["L_end"]),
            " m (floor 1e-6 m), t_converted ",
            fmtv(b["t_converted"] - a["t_converted"]),
            " s, L(900) ",
            fmtv(b["L_900"] - a["L_900"]),
            " m, coating left at 900 s ",
            fmtv(b["coating_frac_900"] - a["coating_frac_900"]),
            ".",
        )
    end

    # ---- grading ------------------------------------------------------------
    println(
        io,
        "\n## Grading: deviation from the finest successful uniform run of each variant (900 s)\n",
    )
    println(
        io,
        "| variant | mesh | stretch | n | ΔTs(100) K | ΔTs(300) K | ΔTs(900) K | ΔTb(900) K | ΔL_end mm | Δpeak MLR g/m²/s | dz0 min µm | dz ratio |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
    grade_qois = (:Ts_100, :Ts_300, :Tb_900, :peak_mlr)
    grading = Dict{String,Vector{NamedTuple}}()
    for v in VARIANTS
        ur = filter(ok, uniform_rows(v))
        isempty(ur) && continue
        ref = ur[end]
        for r in vcat(ur[1:(end-1)], filter(ok, stretched_rows(v)))
            println(
                io,
                @sprintf(
                    "| %s | %s | %.2f | %d | %s | %s | %s | %s | %s | %s | %s | %s |",
                    r["variant"],
                    r["mesh"],
                    r["stretch"],
                    r["n"],
                    fmtv(r["Ts_100"] - ref["Ts_100"]),
                    fmtv(r["Ts_300"] - ref["Ts_300"]),
                    fmtv(r["Ts_900"] - ref["Ts_900"]),
                    fmtv(r["Tb_900"] - ref["Tb_900"]),
                    fmtv(1e3 * (r["L_end"] - ref["L_end"])),
                    fmtv(r["peak_mlr"] - ref["peak_mlr"]),
                    fmtv(1e6 * r["dz0_min"]),
                    fmtv(r["dz_ratio"])
                )
            )
        end
        # Per stretched run: among the QoIs whose uniform error at the same n
        # exceeds the floor, how many does the stretched mesh bring closer?
        for s in filter(ok, stretched_rows(v))
            u = filter(r -> r["n"] == s["n"] && r !== ref, ur)
            isempty(u) && continue
            u = u[1]
            tested =
                [q for q in grade_qois if abs(u[string(q)] - ref[string(q)]) > FLOOR[q]]
            better = count(
                q ->
                    abs(s[string(q)] - ref[string(q)]) < abs(u[string(q)] - ref[string(q)]),
                tested,
            )
            push!(
                get!(grading, v.label, NamedTuple[]),
                (; stretch = s["stretch"], n = s["n"], better, tested = length(tested)),
            )
        end
    end

    # ---- failure probe ------------------------------------------------------
    println(io, "\n## Stretch 1.10, n = 100: where the lateral law throws\n")
    println(
        io,
        "Probe with two process-local recording overrides (copies of `residual!` and of the WithChi ",
        "`update_cross_section_area!`, recording only): the time and scalar type of every residual call and the ",
        "column-averaged χ̄ handed to the lateral law, until the package's positivity check throws. Trace passes: ",
        "`r` = Float64 residual evaluation (a stage predictor or a Newton trial state), `J` = ForwardDiff dual pass ",
        "through `residual!`, `A` = an area update outside `residual!` (the structured Jacobian assembly, ",
        "`src/Jacobian/structured_backend.jl`, calls `update_cross_section_area!` directly at the step's base state; ",
        "it is recognised as a second update tagged with the same residual-call index). ",
        "`accepted` counts the guard callback's fires, i.e. accepted steps; t/dt₀ is the time of the throw over ",
        "`dt_initial` = 1e-3 s, to be read against KenCarp4's stage abscissae c = (0, 1/2, 83/250, 31/50, 17/20, 1). ",
        "Initialization is checked as in the review probe: area, χ range and residual finiteness of the initial state.\n",
    )
    println(
        io,
        "| variant | init area | init χ min/max | init residual finite | outcome | accepted steps | residual calls (assemblies / dual passes / at t = 0) | throw: call, t s, t/dt₀, pass | χ̄ at throw | χ̄ in message | first χ̄ ≠ 0: call, t s, χ̄ | max abs χ̄ before throw | warn | wall s |",
    )
    println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in frows
        println(
            io,
            @sprintf(
                "| %s | %s | %s / %s | %s | %s | %d | %d (%d / %d / %d) | c%d, %s, %s, %s | %s | %s | c%s, %s, %s | %s | %s/%s/%s | %s |",
                r["variant"],
                fmtv(r["init_area"]),
                fmtv(r["init_chi_min"]),
                fmtv(r["init_chi_max"]),
                r["init_residual_finite"],
                first(r["retcode"], 60),
                r["accepted_seen"],
                r["residual_calls"],
                r["jacobian_assemblies"],
                r["dual_passes"],
                r["calls_at_t0"],
                r["residual_calls"],
                fmtv(r["t_throw"]),
                fmtv(r["t_throw_over_dt0"]),
                r["throw_pass"],
                fmtv(r["chi_bar_throw"]),
                fmtv(r["chi_bar_message"]),
                string(r["first_nonzero_call"]),
                fmtv(r["t_first_nonzero"]),
                fmtv(r["chi_bar_first_nonzero"]),
                fmtv(r["chi_bar_max_abs_before_throw"]),
                r["warn_pinned"],
                r["warn_newton"],
                r["warn_other"],
                fmtv(r["wall_s"])
            )
        )
    end
    for r in frows
        trace = split(r["chi_trace"], ";")
        shown =
            length(trace) <= 24 ? trace :
            vcat(trace[1:16], ["…($(length(trace) - 22) calls)…"], trace[(end-5):end])
        println(io, "\n$(r["variant"]) χ̄ trace (call index, t, pass, χ̄):\n")
        println(io, "```\n", join(shown, "\n"), "\n```")
    end

    # ---- findings -----------------------------------------------------------
    println(io, "\n## Findings\n")
    k = 0
    for v in VARIANTS
        c = conv[v.label]
        k += 1
        if isempty(c)
            println(io, "$k. $(v.label): fewer than two successful uniform runs.")
            continue
        end
        ns = c[:Ts_100].ns
        done = [q for q in QOIS if c[q].c !== nothing]
        open_ = [q for q in QOIS if c[q].c === nothing]
        s = "$(v.label) (uniform n = $(join(ns, "/")), finest pair $(ns[end-1])→$(ns[end])): "
        s *=
            isempty(done) ? "no QoI is at its floor" :
            "at the floor: " * join(
                (
                    @sprintf(
                        "%s from n = %d (last Δ = %s %s, floor %s)",
                        q,
                        c[q].c,
                        fmtv(c[q].d[end]),
                        UNIT[q],
                        FLOOR[q]
                    ) for q in done
                ),
                ", ",
            )
        if !isempty(open_)
            s *=
                "; still moving on the finest pair: " * join(
                    (
                        @sprintf(
                            "%s (Δ = %s %s against floor %s, p = %s)",
                            q,
                            fmtv(c[q].d[end]),
                            UNIT[q],
                            FLOOR[q],
                            isempty(c[q].p) ? "–" : fmt_or_dash(c[q].p[end])
                        ) for q in open_
                    ),
                    ", ",
                )
        end
        println(io, "$k. $s.")
    end
    # long horizon
    k += 1
    if length(long_ok) == 2
        a, b = long_ok
        println(
            io,
            "$k. 50× to full conversion: at 900 s the remaining coating mass is ",
            fmtv(a["coating_frac_900"]),
            " (n = $(a["n"])) and ",
            fmtv(b["coating_frac_900"]),
            " (n = $(b["n"])) of the initial, so the 900 s thickness is not a final thickness. ",
            a["converted"] && b["converted"] ?
            @sprintf(
                "Conversion below %s is reached at t = %s s (n = %d) and %s s (n = %d), %s s apart (the stop test runs only on accepted steps, which average %s s over these runs); the thickness at the last saved output before the stop (t = %s s for both) is %s and %s mm, a difference of %s m against the 1e-6 m floor (at 900 s the same pair differs by %s m).",
                CONVERTED_FRAC,
                fmtv(a["t_converted"]),
                a["n"],
                fmtv(b["t_converted"]),
                b["n"],
                fmtv(b["t_converted"] - a["t_converted"]),
                fmtv(
                    round(
                        mean([
                            a["t_converted"] / a["steps"],
                            b["t_converted"] / b["steps"],
                        ]);
                        digits = 1,
                    ),
                ),
                fmtv(a["t_reached"]),
                fmtv(1e3 * a["L_end"]),
                fmtv(1e3 * b["L_end"]),
                fmtv(b["L_end"] - a["L_end"]),
                fmtv(b["L_900"] - a["L_900"])
            ) :
            @sprintf(
                "The horizon of %d s ended before conversion below %s (coating left %s at n = %d, %s at n = %d); the thickness at the end is %s and %s mm (difference %s m).",
                Int(T_LONG),
                CONVERTED_FRAC,
                fmtv(a["coating_frac_end"]),
                a["n"],
                fmtv(b["coating_frac_end"]),
                b["n"],
                fmtv(1e3 * a["L_end"]),
                fmtv(1e3 * b["L_end"]),
                fmtv(b["L_end"] - a["L_end"])
            ),
        )
    else
        println(
            io,
            "$k. 50× to full conversion: ",
            length(long_ok),
            " of ",
            length(lr),
            " long-horizon runs finished.",
        )
    end
    # grading
    k += 1
    gparts = String[]
    for v in VARIANTS
        g = get(grading, v.label, NamedTuple[])
        isempty(g) && continue
        push!(
            gparts,
            "$(v.label): " * join(
                (
                    @sprintf(
                        "stretch %.2f n = %d closer on %d of %d",
                        x.stretch,
                        x.n,
                        x.better,
                        x.tested
                    ) for x in g
                ),
                ", ",
            ),
        )
    end
    println(
        io,
        "$k. Stretched versus uniform at the same n, counted over the QoIs among Ts(100), Ts(300), Tb(900) and peak MLR whose uniform error against the finest uniform run exceeds the floor: ",
        isempty(gparts) ? "no comparison possible" : join(gparts, "; "),
        ".",
    )
    # failures
    k += 1
    bad = filter(r -> !ok(r) && !(r["retcode"] == "Terminated" && r["converted"]), rows)
    if isempty(bad)
        println(
            io,
            "$k. All $(length(rows)) runs returned `Success` (or stopped at conversion) within the $(Int(TIMEOUT_S)) s guard.",
        )
    else
        println(
            io,
            "$k. ",
            length(bad),
            " of ",
            length(rows),
            " runs did not succeed: ",
            join(
                (
                    @sprintf(
                        "%s %s stretch=%.2f n=%d (smallest initial cell %s µm, cell-size ratio %s) → %s%s; accepted steps %d, last accepted t = %s s, wall %s s",
                        r["variant"],
                        r["mesh"],
                        r["stretch"],
                        r["n"],
                        fmtv(1e6 * r["dz0_min"]),
                        fmtv(r["dz0_max"] / r["dz0_min"]),
                        first(r["retcode"], 110),
                        r["timed_out"] ? " (timed out)" : "",
                        r["accepted_seen"],
                        fmtv(r["t_last_accepted"]),
                        fmtv(r["wall_s"])
                    ) for r in bad
                ),
                "; ",
            ),
            ". Successful runs took ",
            fmtv(minimum(r["wall_s"] for r in rows if ok(r))),
            "–",
            fmtv(maximum(r["wall_s"] for r in rows if ok(r))),
            " s each.",
        )
    end
    k += 1
    if !isempty(frows)
        ok_init = all(
            r ->
                r["init_area"] == 1.0 &&
                r["init_chi_max"] == 0.0 &&
                r["init_residual_finite"],
            frows,
        )
        no_steps = all(r -> r["accepted_seen"] == 0, frows)
        println(
            io,
            "$k. Re-diagnosis of the stretch-1.10, 100-cell failures (smallest initial cell 22 nm, cell-size ratio 1.25e4): ",
            ok_init ?
            "all three decks initialize with cross-section area 1.0, χ = 0 in every cell and a finite residual, so this is not an initialization failure. " :
            "initialization differs between the decks (see the table). ",
            no_steps ? "No step is accepted before the throw. " :
            "At least one step is accepted before the throw. ",
            join(
                (
                    @sprintf(
                        "%s: the law throws on residual call %d at t = %s s = %s dt₀ (pass %s), χ̄ = %s, after %d residual call(s) at t = 0, %d structured Jacobian assemblies and %d dual passes; the first non-zero χ̄ appears on call %s (t = %s s, χ̄ = %s) and the largest |χ̄| before the throw is %s",
                        r["variant"],
                        r["residual_calls"],
                        fmtv(r["t_throw"]),
                        fmtv(r["t_throw_over_dt0"]),
                        r["throw_pass"],
                        fmtv(r["chi_bar_throw"]),
                        r["calls_at_t0"],
                        r["jacobian_assemblies"],
                        r["dual_passes"],
                        string(r["first_nonzero_call"]),
                        fmtv(r["t_first_nonzero"]),
                        fmtv(r["chi_bar_first_nonzero"]),
                        fmtv(r["chi_bar_max_abs_before_throw"])
                    ) for r in frows
                ),
                "; ",
            ),
            ". The times of the residual calls are 0, 0.5 dt₀, 0.332 dt₀ and 0.62 dt₀: the stage abscissae c₂, c₃, c₄ of the first KenCarp4 step at `dt_initial`, so the χ̄ that reaches the lateral law is a Newton trial state of the first step's stage at c₄ = 31/50, not the initial state. ",
            "The stages at c₂ and c₃ converge with χ̄ of order 1e-18 to 1e-6 (the trace is identical for the three expansion factors up to the throw, since no conversion has happened yet); at c₄ the first Newton correction (call 13) takes χ̄ to 8.4e83 and the second (call 14) to a negative value of order 1e68 to 1e70, where 1 + 0.3 χ̄ < 0 and the package's positivity check throws. ",
            "This is a first-step trial-state failure of the lateral law (the law is evaluated on an unbounded Newton iterate, with no clamp of χ̄ to [0, 1] and no step rejection before the throw), not an initialization failure. The three surface-temperature clamp and two surface-Newton warnings belong to the same trial states.",
        )
    end
end
println("wrote ", CSV, ", ", FAIL_CSV, " and ", MD)
println(@sprintf("total script wall time %.1f s", time() - t_script))
