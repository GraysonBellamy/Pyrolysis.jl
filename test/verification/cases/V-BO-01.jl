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
#
# On the finite-volume core the mesh moves by the equidistribution correction,
# the PMMA–Kaowool interface a Lagrangian node, so the PMMA cells thin
# together and the layer deactivates in one invocation into the Kaowool cell
# beneath; the surface QoI there counts accepted states whose surface balance
# did not converge or sits on the trial clamp, since that Newton never warns.

# The deck, the warning counters and the step recorder below are shared by
# V-BO-02, V-BO-03 and V-BO-05 (their files load after this one).

using Logging
using OrdinaryDiffEq: CallbackSet, DiscreteCallback, SciMLBase

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
    A = 1.35e11,
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
    rx = Reaction(:pmma_gasification, 1 => 3; A, E = 164e3, h = 8.7e5, n = 1.0)
    return Material(
        name = :PMMA_MaCFP,
        components = (pmma, kw, gas),
        reactions = (rx,),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
end

# `n_pmma` cells across the PMMA layer; the Kaowool takes the same pitch.
# Surface-absorption runs pass zero absorption coefficients (unused there);
# the gradient gate perturbs the absorption coefficient, the incident flux
# and the pre-exponential factor `A`.
function _vbo_problem(
    n_pmma::Int;
    substrate_alpha = 8.0e6 / VBO_RHO_KW,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
    flux = 50e3,
    A = 1.35e11,
    t_end = 500.0,
)
    mat = _vbo_material(; substrate_alpha, pmma_alpha, A)
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
        RadiativeFluxBC(; flux, absorptivity = 0.96) +
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

# What the step recorder saw over the accepted steps of a solve: the smallest
# step, the largest live-cell temperature, the time of the last step (a stop
# between output times is not saved), on the new core the accepted states
# whose surface energy balance did not converge or sits on the trial clamp of
# the property bracket, and the deactivation threshold (`NaN` on the old core,
# whose surface Newton warns instead).
mutable struct _VBORecord
    dt_min::Float64
    T_max::Float64
    t_last::Float64
    clamped::Int
    threshold::Float64
end
_VBORecord() = _VBORecord(Inf, -Inf, NaN, 0, NaN)

# Warning counters and step records of the most recent solve of each case
# (the QoI closures read them after `solve_fn` returns).
const VBO_SURFACE_WARNINGS = Dict{String,Int}()
const VBO_TRANSMITTED_WARNINGS = Dict{String,Int}()
const VBO_VOID_WARNINGS = Dict{String,Int}()
const VBO_STEP_STATS = Dict{String,_VBORecord}()

# The step recorder: after every accepted step, the step and its time, then
# `observe!(u, t)` on the accepted state.
function _vbo_recorder(observe!, rec::_VBORecord)
    return DiscreteCallback(
        (u, t, integ) -> true,
        integ -> begin
            rec.dt_min = min(rec.dt_min, integ.t - integ.tprev)
            observe!(integ.u, integ.t)
            rec.t_last = integ.t
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
end

# The live boundary faces of the thermal kernel whose surface energy balance,
# at the last residual evaluation, did not converge or ended where the trial
# clamp departs from the identity. The new core's surface Newton evaluates the
# condition at the clamped iterate and never warns: on a trial state the clamp
# is its purpose, on an accepted state it would be a failure.
function _vbo_clamped(thermal, ds, topo)
    n = 0
    for p in thermal.patches, f0 in p.faces
        f = Pyrolysis.Meshes.live_face(topo, ds.active, f0)
        f == 0 && continue
        s = thermal.states[f]
        clamped = Pyrolysis.States.trial_temperature(s.T_surface) != s.T_surface
        (s.converged && !clamped) || (n += 1)
    end
    return n
end

# The PMMA–Kaowool interface node of a V-BO deck, from the back face: the
# Kaowool takes the cells of its thickness at the mesh's uniform pitch.
function _vbo_interface_node(problem)
    nodes = problem.mesh.nodes
    dz = (nodes[end].z - nodes[1].z) / length(problem.mesh.cells)
    return round(Int, VBO_L_KW / dz) + 1
end

# The new core's solve with the step recorder, which evaluates the residual
# at every accepted state and reads the live cells' temperatures and the
# surface balances. The old core's depletion keywords have no counterpart
# there. `instrument(cp)` returns a further callback of the set-up problem
# for the gate scripts (`nothing` for none); `min_thickness` and, when given,
# `located_events` go to the integration, the rest to the setup.
function _vbo_core_solve(
    problem,
    rec::_VBORecord;
    handle_depletion = false,
    depletion_threshold = nothing,
    min_cells = nothing,
    min_thickness = 1e-6,
    located_events = nothing,
    instrument = cp -> nothing,
    kwargs...,
)
    cp = Pyrolysis.Integration.setup(
        problem;
        energy_state = VERIFICATION_ENERGY_STATE,
        deactivation = handle_depletion,
        kwargs...,
    )
    def, ws = cp.def, cp.ws
    topo = def.setup.topology
    thermal =
        Pyrolysis.Kernels.find_kernel(Pyrolysis.Kernels.ThermalBoundaryKernel, def.kernels)
    du = similar(cp.u0)
    observe!(u, t) = begin
        Pyrolysis.Assembly.residual!(du, u, def, t, ws)
        ds, _ = Pyrolysis.Assembly.buffers(ws, def, Float64)
        rec.T_max =
            max(rec.T_max, maximum(ds.T[i] for i in eachindex(ds.T) if ds.active[i]))
        rec.clamped += _vbo_clamped(thermal, ds, topo)
    end
    cp.deactivation === nothing || (rec.threshold = cp.deactivation.threshold)
    extra = instrument(cp)
    recorder = _vbo_recorder(observe!, rec)
    callback = extra === nothing ? recorder : CallbackSet(recorder, extra)
    located = located_events === nothing ? (;) : (; located_events)
    return Pyrolysis.Integration.integrate(cp; callback, min_thickness, located...)
end

# The case's solve. On the new core a `corrected` case moves its mesh by the
# equidistribution correction, with the layer interface Lagrangian (given
# keywords take precedence).
function _vbo_solve_fn(id::String; corrected::Bool = false)
    return (problem; kwargs...) -> begin
        surf = _VBOWarnCounter(Logging.current_logger(), r"^Surface", Ref(0))
        tr = _VBOWarnCounter(surf, r"transmitted through the back face", Ref(0))
        void = _VBOWarnCounter(tr, r"gasified from within", Ref(0))
        # On the old core the maximum runs over every temperature slot: a
        # ghost keeps the temperature it was parked with, which lies between
        # two live-cell temperatures.
        n = length(problem.mesh.cells)
        rec = _VBORecord()
        sol = Logging.with_logger(void) do
            if VERIFICATION_CORE == "new"
                kw =
                    corrected ?
                    (;
                        correction = Pyrolysis.Closures.EquidistributionCorrection(),
                        lagrangian = [_vbo_interface_node(problem)],
                        kwargs...,
                    ) : kwargs
                _vbo_core_solve(problem, rec; kw...)
            else
                recorder = _vbo_recorder(
                    (u, t) -> (rec.T_max = max(rec.T_max, maximum(view(u, 1:n)))),
                    rec,
                )
                solve(problem; progress = false, callback = recorder, kwargs...)
            end
        end
        VBO_SURFACE_WARNINGS[id] = surf.count[]
        VBO_TRANSMITTED_WARNINGS[id] = tr.count[]
        VBO_VOID_WARNINGS[id] = void.count[]
        VBO_STEP_STATS[id] = rec
        return sol
    end
end

# The peak and its time are the main peak's: over the saved times before the
# solve's last deactivation event, whose merge (or park) releases the layer's
# residue in a burst the gate judges by refinement instead.
function _vbo_measures(sol)
    mlr = sol.extras.gas_generation_rate .* 1000           # g m⁻² s⁻¹
    parks = sol.extras.park_times
    before =
        isempty(parks) ? length(sol.t) : something(findlast(<(maximum(parks)), sol.t), 1)
    ipk = argmax(view(mlr, 1:before))
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
    ledger_tol = 1e-3,
    ledger_comment = "total-mass drop vs ∫ ṁ″ dt over the run",
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
            tol = ledger_tol,
            comment = ledger_comment,
        ),
    )
    surface_warnings && push!(out, _vbo_surface_qoi(id))
    return out
end

# The surface energy balance's QoI: on the old core no warning from its
# Newton (converged, never pinned at the clamp) at any evaluation; on the new
# core, whose Newton clamps its trial iterates smoothly and never warns, no
# accepted state whose balance did not converge or sits on the clamp.
function _vbo_surface_qoi(id)
    if VERIFICATION_CORE == "new"
        n = get(() -> _VBORecord(-1.0, -1.0, NaN, -1, NaN), VBO_STEP_STATS, id).clamped
        return QoIResult(
            name = "clamped surface balances on accepted states",
            value = Float64(n),
            exact = 0.0,
            error = Float64(abs(n)),
            tol = 0.5,
            comment = "boundary temperature solve converged inside the property bracket at every accepted step",
        )
    end
    n_warn = get(VBO_SURFACE_WARNINGS, id, -1)
    return QoIResult(
        name = "surface-Newton warnings",
        value = Float64(n_warn),
        exact = 0.0,
        error = Float64(n_warn),
        tol = 0.5,
        comment = "boundary temperature solve converged at every evaluation",
    )
end

# The deactivation invocations of a solve, each with its cells below the
# threshold `ε_m`, the cells the band took and its largest merged fraction.
function _vbo_invocations(sol, ε_m)
    return [
        begin
            es = filter(e -> e.invocation == k, sol.events)
            below = count(e -> e.fraction < ε_m, es)
            (;
                invocation = k,
                t = first(es).t,
                cells = length(es),
                below,
                band = length(es) - below,
                largest = maximum(e.fraction for e in es),
            )
        end for k in unique(e.invocation for e in sol.events)
    ]
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
        VERIFICATION_CORE == "new" || return out
        # The correction carries the layer to one deactivation, and the band
        # merges no cell that still holds more than 5 ε_m of its mass.
        ε_m = VBO_STEP_STATS[id].threshold
        inv = _vbo_invocations(sol, ε_m)
        largest = maximum(v.largest for v in inv; init = 0.0)
        push!(
            out,
            QoIResult(
                name = "one invocation for the PMMA layer",
                value = Float64(length(inv)),
                exact = 1.0,
                error = abs(length(inv) - 1.0) + (largest <= 5ε_m ? 0.0 : 1.0),
                tol = 0.5,
                comment = "$(sum(v.cells for v in inv; init = 0)) cells; largest merged fraction $(round(largest; sigdigits = 3)) against 5 ε_m = $(5ε_m)",
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
        solve_fn = _vbo_solve_fn(id; corrected = true),
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
