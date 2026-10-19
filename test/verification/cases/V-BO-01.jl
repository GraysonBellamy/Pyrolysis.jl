# =============================================================================
# V-BO-01 — PMMA-on-Kaowool burnout under ALE (MaCFP gasification deck)
# =============================================================================
#
# A 5.9 mm PMMA slab on 12 mm Kaowool under 50 kW/m² burns out completely:
# every PMMA cell shrinks to the depletion threshold and is parked, the
# exposed boundary walks down the stack, and the run ends on the bare
# substrate. The deck is the MaCFP PMMA gasification case with the Gpyro V2
# property fits (power-law k(T), c(T) above 394 K — positive on the whole
# property bracket) and in-depth (Beer–Lambert) absorption.
#
# The run must burn out completely, lose exactly the mass of the PMMA layer,
# match the reference peak mass loss rate, its time and the burnout time,
# close the mass ledger, and converge the surface temperature at every
# evaluation, within a wall-time budget at 0.098 mm cells.

# The deck, the warning counters and the step recorder below are shared by
# V-BO-02, V-BO-03 and V-BO-05 (their files load after this one).

using Logging
using OrdinaryDiffEq: DiscreteCallback, SciMLBase

const VBO_L_PMMA = 0.0059
const VBO_L_KW = 0.012
const VBO_T0 = 300.0
const VBO_T_INF = 298.0
const VBO_RHO_PMMA = 1210.0
const VBO_RHO_KW = 256.0
const VBO_MASS_LOSS_EXACT =
    100 * VBO_RHO_PMMA * VBO_L_PMMA / (VBO_RHO_PMMA * VBO_L_PMMA + VBO_RHO_KW * VBO_L_KW)

function _vbo_material(;
    substrate_alpha = 8.0e6 / VBO_RHO_KW,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
)
    k_pmma = CallableProperty(T -> T <= 394.0 ? 0.16 : 0.2156 * (T / 300.0)^-0.745)
    c_pmma = CallableProperty(
        T ->
            T <= 394.0 ? 1127.0793 * (T / 300.0)^1.9283 : 1768.9067 * (T / 300.0)^0.5552,
    )
    pmma = SolidComponent(
        :PMMA,
        ρ = VBO_RHO_PMMA,
        c = c_pmma,
        k = k_pmma,
        ε = 0.96,
        α = pmma_alpha,
        λ = 2e-5,
    )
    kw = SolidComponent(
        :KAOWOOL,
        ρ = VBO_RHO_KW,
        c = 1100.0,
        k = 0.058,
        ε = 1.0,
        λ = 2e-5,
        α = substrate_alpha,
    )
    gas = GaseousComponent(:MMA, M = 0.044, c = 0.0, k = 0.02, λ = 2e-2)
    rx = Reaction(:pmma_gasification, 1 => 3; A = 1.35e11, E = 164e3, h = 8.7e5, n = 1.0)
    return Material(
        name = :PMMA_MaCFP,
        components = (pmma, kw, gas),
        reactions = (rx,),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
end

# `n_pmma` cells across the PMMA layer; the Kaowool takes the same pitch.
# Surface-absorption runs pass zero absorption coefficients (unused there).
function _vbo_problem(
    n_pmma::Int;
    substrate_alpha = 8.0e6 / VBO_RHO_KW,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
    t_end = 500.0,
)
    mat = _vbo_material(; substrate_alpha, pmma_alpha)
    dz = VBO_L_PMMA / n_pmma
    n_kw = round(Int, VBO_L_KW / dz)
    n = n_pmma + n_kw
    mesh = create_uniform_mesh(
        n * dz,
        n,
        Val(3);
        T_initial = VBO_T0,
        ξ_initial = (VBO_RHO_PMMA, 0.0, 0.0),
    )
    z_if = n_kw * dz
    ξ0 =
        [z -> (z > z_if ? VBO_RHO_PMMA : 0.0), z -> (z > z_if ? 0.0 : VBO_RHO_KW), z -> 0.0]
    top =
        RadiativeFluxBC(flux = 50e3, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = VBO_T_INF) +
        ConvectiveBC(h = 8.2, T_inf = VBO_T_INF)
    bottom = RadiativeBC(ε = 1.0, T_inf = VBO_T_INF)
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = bottom),
        mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bcs,
        T_initial = z -> VBO_T0,
        ξ_initial = ξ0,
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:0.5:t_end),
    )
end

# Counts warnings whose message matches `pattern` while forwarding everything
# to the enclosing logger. The surface-Newton QoIs read the count.
struct _VBOWarnCounter{L<:Logging.AbstractLogger} <: Logging.AbstractLogger
    inner::L
    pattern::Regex
    count::Base.RefValue{Int}
end
Logging.min_enabled_level(l::_VBOWarnCounter) = Logging.min_enabled_level(l.inner)
Logging.shouldlog(l::_VBOWarnCounter, args...) = Logging.shouldlog(l.inner, args...)
Logging.catch_exceptions(l::_VBOWarnCounter) = Logging.catch_exceptions(l.inner)
function Logging.handle_message(l::_VBOWarnCounter, level, message, args...; kwargs...)
    if level >= Logging.Warn && occursin(l.pattern, string(message))
        l.count[] += 1
    end
    return Logging.handle_message(l.inner, level, message, args...; kwargs...)
end

# Warning counters and step statistics of the most recent solve of each case
# (the QoI closures read them after `solve_fn` returns).
const VBO_SURFACE_WARNINGS = Dict{String,Int}()
const VBO_TRANSMITTED_WARNINGS = Dict{String,Int}()
const VBO_VOID_WARNINGS = Dict{String,Int}()
const VBO_STEP_STATS =
    Dict{String,@NamedTuple{dt_min::Float64,T_max::Float64,t_last::Float64}}()

# The step recorder: the smallest accepted step, the largest temperature and
# the time of the last accepted step (a stop between output times is not
# saved), with the temperatures of a state read by `temperatures(u)`.
function _vbo_recorder(temperatures, dt_min, T_max, t_last)
    return DiscreteCallback(
        (u, t, integ) -> true,
        integ -> begin
            dt_min[] = min(dt_min[], integ.t - integ.tprev)
            T_max[] = max(T_max[], temperatures(integ.u))
            t_last[] = integ.t
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
end

# The new core's solve with the step recorder, which reads the live cells'
# temperatures through the layout and the state mapping. The old core's
# depletion keywords have no counterpart there.
function _vbo_core_solve(
    problem,
    dt_min,
    T_max,
    t_last;
    handle_depletion = false,
    depletion_threshold = nothing,
    min_cells = nothing,
    min_thickness = 1e-6,
    kwargs...,
)
    cp = Pyrolysis.Integration.setup(
        problem;
        energy_state = VERIFICATION_ENERGY_STATE,
        deactivation = handle_depletion,
        kwargs...,
    )
    layout, setup, ds = cp.def.layout, cp.def.setup, cp.ws.derived
    live_max(u) = begin
        Pyrolysis.States.derive!(ds, u, layout, setup)
        maximum(ds.T[i] for i in eachindex(ds.T) if ds.active[i])
    end
    recorder = _vbo_recorder(live_max, dt_min, T_max, t_last)
    return Pyrolysis.Integration.integrate(cp; callback = recorder, min_thickness)
end

function _vbo_solve_fn(id::String)
    return (problem; kwargs...) -> begin
        surf = _VBOWarnCounter(Logging.current_logger(), r"^Surface", Ref(0))
        tr = _VBOWarnCounter(surf, r"transmitted through the back face", Ref(0))
        void = _VBOWarnCounter(tr, r"gasified from within", Ref(0))
        # On the old core the maximum runs over every temperature slot: a
        # ghost keeps the temperature it was parked with, which lies between
        # two live-cell temperatures.
        n = length(problem.mesh.cells)
        dt_min = Ref(Inf)
        T_max = Ref(-Inf)
        t_last = Ref(NaN)
        sol = Logging.with_logger(void) do
            if VERIFICATION_CORE == "new"
                _vbo_core_solve(problem, dt_min, T_max, t_last; kwargs...)
            else
                recorder = _vbo_recorder(u -> maximum(view(u, 1:n)), dt_min, T_max, t_last)
                solve(problem; progress = false, callback = recorder, kwargs...)
            end
        end
        VBO_SURFACE_WARNINGS[id] = surf.count[]
        VBO_TRANSMITTED_WARNINGS[id] = tr.count[]
        VBO_VOID_WARNINGS[id] = void.count[]
        VBO_STEP_STATS[id] = (; dt_min = dt_min[], T_max = T_max[], t_last = t_last[])
        return sol
    end
end

function _vbo_measures(sol)
    mlr = sol.extras.gas_generation_rate .* 1000           # g m⁻² s⁻¹
    ipk = argmax(mlr)
    th = sol.extras.thickness_history
    i_bo = findfirst(<(VBO_L_KW + 0.02e-3), th)
    t_bo = i_bo === nothing ? NaN : sol.t[i_bo]
    m = sol.extras.total_mass
    mass_loss = 100 * (1 - m[end] / m[1])
    # Mass ledger: total-mass drop vs ∫ ṁ″ dt (trapezoid on output times).
    ṁ = sol.mass_loss_rate
    gen = sum(0.5 * (ṁ[j] + ṁ[j+1]) * (sol.t[j+1] - sol.t[j]) for j = 1:(length(sol.t)-1))
    ledger = abs((m[1] - m[end]) - gen) / (m[1] - m[end])
    return (;
        mlr,
        peak = mlr[ipk],
        t_peak = sol.t[ipk],
        t_bo,
        mass_loss,
        ledger,
        L_end = th[end],
    )
end

# Maximum live-cell temperature at each output time (ghost slots read NaN).
_vbo_max_live_T(sol) =
    [maximum(x -> isnan(x) ? -Inf : x, view(sol.T, :, k)) for k in eachindex(sol.t)]

# Report figure of one time history.
_vbo_history_curve(id, name, title, ylabel, label, t, y) = ReportCurve(
    name = name,
    title = "$(id): $(title)",
    xlabel = "t (s)",
    ylabel = ylabel,
    series = [CurveSeries(label = label, x = collect(Float64, t), y = collect(Float64, y))],
)

# Peak and burnout time of the corrected-deck reference run (0.098 mm cells).
const VBO_PEAK_REF = 32.39      # g m⁻² s⁻¹
const VBO_T_PEAK_REF = 335.0    # s
const VBO_T_BO_REF = 363.5      # s (live thickness < 12.02 mm)

function _vbo_common_qois(
    m;
    id,
    peak_tol,
    peak_ref = VBO_PEAK_REF,
    surface_warnings = true,
    peak_comment = "g m⁻² s⁻¹, gas generation rate",
)
    out = QoIResult[]
    push!(
        out,
        QoIResult(
            name = "complete burnout: live thickness at t_end is the substrate",
            value = m.L_end,
            exact = VBO_L_KW,
            error = abs(m.L_end - VBO_L_KW),
            tol = 1e-5,
            comment = "every PMMA cell parked; 12.00 ± 0.01 mm",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "mass loss of the gasified PMMA layer",
            value = m.mass_loss,
            exact = VBO_MASS_LOSS_EXACT,
            error = abs(m.mass_loss - VBO_MASS_LOSS_EXACT),
            tol = 0.3,
            comment = "percent of the initial mass; exact = ρ_PMMA L_PMMA / m₀",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "peak MLR",
            value = m.peak,
            exact = peak_ref,
            error = abs(m.peak - peak_ref) / peak_ref,
            tol = peak_tol,
            comment = peak_comment,
        ),
    )
    push!(
        out,
        QoIResult(
            name = "mass ledger closure",
            value = m.ledger,
            exact = 0.0,
            error = m.ledger,
            tol = 1e-3,
            comment = "total-mass drop vs ∫ ṁ″ dt over the run",
        ),
    )
    if surface_warnings
        n_warn = get(VBO_SURFACE_WARNINGS, id, -1)
        push!(
            out,
            QoIResult(
                name = "surface-Newton warnings",
                value = Float64(n_warn),
                exact = 0.0,
                error = Float64(n_warn),
                tol = 0.5,
                comment = "boundary temperature solve converged at every evaluation",
            ),
        )
    end
    return out
end

function _verification_case_vbo01()
    id = "V-BO-01"
    function qois(sol)
        m = _vbo_measures(sol)
        out = _vbo_common_qois(m; id, peak_tol = 0.03)
        push!(
            out,
            QoIResult(
                name = "time of peak MLR",
                value = m.t_peak,
                exact = VBO_T_PEAK_REF,
                error = abs(m.t_peak - VBO_T_PEAK_REF) / VBO_T_PEAK_REF,
                tol = 0.03,
                comment = "s",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "burnout time (live thickness < 12.02 mm)",
                value = m.t_bo,
                exact = VBO_T_BO_REF,
                error = abs(m.t_bo - VBO_T_BO_REF) / VBO_T_BO_REF,
                tol = 0.01,
                comment = "s",
            ),
        )
        return out
    end

    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout under ALE (parked depletion)",
        tags = [:ale, :burnout, :depletion, :beer_lambert, :coupled, :performance],
        references = [
            "MaCFP-2 PMMA gasification case (NIST), 50 kW/m²",
            "Gpyro V2 PMMA deck property fits",
        ],
        build = n -> _vbo_problem(n),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = (
            use_ale = true,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
            min_thickness = 1e-3,
            radiation_model = BEER_LAMBERT,
            abstol = 1e-8,
            reltol = 1e-6,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 60, budget_s = 20.0),
        report = (
            n_cells_ladder = [30, 60, 120],
            qoi_curves = (:mlr_history, :thickness_history),
            curves = sol -> begin
                m = _vbo_measures(sol)
                names, cols = read_reference_csv(
                    "gpyro/vbo01_gpyro_v2_pmma.csv";
                    header_rows = 1,
                )
                t_g = cols[1]
                mlr_g = cols[2]                      # Gpyro writes MLR in g m⁻² s⁻¹
                [
                    ReportCurve(
                        name = :mlr_history,
                        title = "V-BO-01: PMMA gasification rate through burnout",
                        xlabel = "t (s)",
                        ylabel = "MLR (g m⁻² s⁻¹)",
                        series = [
                            CurveSeries(
                                label = "gas generation rate",
                                x = sol.t,
                                y = m.mlr,
                            ),
                            CurveSeries(
                                label = "Gpyro V2",
                                x = t_g,
                                y = mlr_g,
                                crosscode = true,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :thickness_history,
                        title = "V-BO-01: live thickness",
                        xlabel = "t (s)",
                        ylabel = "thickness (m)",
                        series = [
                            CurveSeries(
                                label = "live thickness",
                                x = sol.t,
                                y = sol.extras.thickness_history,
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo01)
