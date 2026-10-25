# =============================================================================
# V-BO-01 — PMMA-on-Kaowool burnout under ALE (MaCFP gasification deck)
# =============================================================================
#
# A 5.9 mm PMMA slab on 12 mm Kaowool under 50 kW/m² burns out completely:
# the PMMA cells thin to the deactivation threshold, the layer is deactivated
# and the run ends on the bare substrate. The deck is the MaCFP PMMA
# gasification case with the Gpyro V2 property fits (power-law k(T), c(T)
# above 394 K — positive on the whole property bracket) and in-depth
# (Beer–Lambert) absorption.
#
# The run must burn out completely, lose exactly the mass of the PMMA layer,
# close the mass ledger, and converge the surface temperature at every
# accepted step, within a wall-time budget at 0.098 mm cells.
# Peak mass loss rate and timing are reported through the history curves;
# the report also compares the histories with Gpyro V2.
#
# The mesh moves by the equidistribution correction, the PMMA–Kaowool
# interface a Lagrangian node, so the PMMA cells thin together and the layer
# deactivates in one invocation into the Kaowool cell beneath. The surface
# QoI counts accepted states whose surface balance did not converge or sits
# on the trial clamp, since that Newton never warns.

# The deck, the step recorder and the QoIs below are shared by V-BO-02,
# V-BO-03 and V-BO-05 to V-BO-08 (their files load after this one).

using OrdinaryDiffEq: CallbackSet, ContinuousCallback, DiscreteCallback, SciMLBase

const VBO_L_PMMA = 0.0059
const VBO_L_KW = 0.012
const VBO_T0 = 300.0
const VBO_T_INF = 298.0
const VBO_RHO_PMMA = 1210.0
const VBO_RHO_KW = 256.0
const VBO_MASS_LOSS_EXACT =
    100 * VBO_RHO_PMMA * VBO_L_PMMA / (VBO_RHO_PMMA * VBO_L_PMMA + VBO_RHO_KW * VBO_L_KW)

# The PMMA of the deck (the Gpyro V2 fits) under `name` and density `ρ`, the
# Kaowool substrate and the MMA gas.
function _vbo_pmma(;
    name = :PMMA,
    ρ = VBO_RHO_PMMA,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
    porosity = 0.01,
)
    k_pmma = CallableProperty(T -> T <= 394.0 ? 0.16 : 0.2156 * (T / 300.0)^-0.745)
    c_pmma = CallableProperty(
        T ->
            T <= 394.0 ? 1127.0793 * (T / 300.0)^1.9283 : 1768.9067 * (T / 300.0)^0.5552,
    )
    return SolidComponent(
        name,
        ρ = ρ,
        φ = porosity,
        c = c_pmma,
        k = k_pmma,
        ε = 0.96,
        α = pmma_alpha,
        λ = 2e-5,
    )
end
_vbo_kaowool(substrate_alpha) = SolidComponent(
    :KAOWOOL,
    ρ = VBO_RHO_KW,
    φ = 0.8,
    c = 1100.0,
    k = 0.058,
    ε = 1.0,
    λ = 2e-5,
    α = substrate_alpha,
)
_vbo_gas() = GaseousComponent(:MMA, M = 0.10012, c = 0.0, k = 0.02)

function _vbo_material(;
    substrate_alpha = 8.0e6 / VBO_RHO_KW,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
    A = 1.35e11,
    porosity = 0.01,
)
    rx = Reaction(:pmma_gasification, 1 => 3; A, E = 164e3, h = 8.7e5, n = 1.0)
    return Material(
        name = :PMMA_MaCFP,
        components = (
            _vbo_pmma(; pmma_alpha, porosity),
            _vbo_kaowool(substrate_alpha),
            _vbo_gas(),
        ),
        reactions = (rx,),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
end

# The node positions of `n_pmma` PMMA cells across `l_pmma` on the Kaowool,
# each PMMA cell `1/r` of the one beneath it, the Kaowool at the pitch of the
# PMMA cell beside it, with the Kaowool's cell count.
function _vbo_nodes(n_pmma, l_pmma, r)
    h = l_pmma * (1 - 1 / r) / (1 - r^(-n_pmma))
    n_kw = round(Int, VBO_L_KW / h)
    pmma = cumsum([h * r^(1 - k) for k = 1:n_pmma])
    return [[k * h for k = 0:n_kw]; n_kw * h .+ pmma], n_kw
end

# The deck's conditions on `mesh` with `material` and initial concentrations
# `ξ0`, saved every `save_every` seconds.
function _vbo_deck(mesh, material, ξ0; flux = 50e3, t_end = 500.0, save_every = 0.5)
    top =
        RadiativeFluxBC(; flux, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = VBO_T_INF) +
        ConvectiveBC(h = 8.2, T_inf = VBO_T_INF)
    bottom = RadiativeBC(ε = 1.0, T_inf = VBO_T_INF)
    bcs = BoundaryConditionSet(
        thermal = (top = top, bottom = bottom),
        mass = (
            top = AmbientGasBC(GasReservoir(P = 0.0, T = VBO_T_INF)),
            bottom = ImpermeableBC(),
        ),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bcs,
        T_initial = z -> VBO_T0,
        ξ_initial = ξ0,
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(0.0:save_every:t_end),
    )
end

# `n_pmma` cells across the PMMA layer of thickness `l_pmma`, uniform or graded
# toward the surface by `grading` (`_vbo_nodes`); the Kaowool takes the pitch
# of the PMMA beside it. Surface-absorption runs pass zero absorption
# coefficients (unused there); the gradient gate perturbs the absorption
# coefficient, the incident flux and the pre-exponential factor `A`.
function _vbo_problem(
    n_pmma::Int;
    substrate_alpha = 8.0e6 / VBO_RHO_KW,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
    flux = 50e3,
    A = 1.35e11,
    t_end = 500.0,
    save_every = 0.5,
    l_pmma = VBO_L_PMMA,
    grading = 1.0,
    porosity = 0.01,
)
    mat = _vbo_material(; substrate_alpha, pmma_alpha, A, porosity)
    if grading == 1
        dz = l_pmma / n_pmma
        n_kw = round(Int, VBO_L_KW / dz)
        n = n_pmma + n_kw
        mesh = create_uniform_mesh(n * dz, n)
        z_if = n_kw * dz
    else
        nodes, n_kw = _vbo_nodes(n_pmma, l_pmma, grading)
        mesh = create_mesh_from_nodes(nodes)
        z_if = nodes[n_kw+1]
    end
    ξ0 =
        [z -> (z > z_if ? VBO_RHO_PMMA : 0.0), z -> (z > z_if ? 0.0 : VBO_RHO_KW), z -> 0.0]
    return _vbo_deck(mesh, mat, ξ0; flux, t_end, save_every)
end

# What the step recorder saw over the accepted steps of a solve: the smallest
# step, the largest live-cell temperature, the time of the last step (a stop
# between output times is not saved), the accepted states whose surface
# energy balance did not converge or sits on the trial clamp of the property
# bracket, the deactivation threshold and the mass ledger.
mutable struct _VBORecord
    dt_min::Float64
    T_max::Float64
    t_last::Float64
    clamped::Int
    threshold::Float64
    ledger::Float64
end
_VBORecord() = _VBORecord(Inf, -Inf, NaN, 0, NaN, NaN)

# Step records of the most recent solve of each case (the QoI closures read
# them after `solve_fn` returns).
const VBO_STEP_STATS = Dict{String,_VBORecord}()

# The step recorder: after every accepted step and its events, the step and
# its time, then `observe!(integ)`.
function _vbo_recorder(observe!, rec::_VBORecord)
    return DiscreteCallback(
        (u, t, integ) -> true,
        integ -> begin
            rec.dt_min = min(rec.dt_min, integ.t - integ.tprev)
            observe!(integ)
            rec.t_last = integ.t
            SciMLBase.derivative_discontinuity!(integ, false)
        end;
        save_positions = (false, false),
    )
end

# The live boundary faces of the thermal kernel whose surface energy balance,
# at the last residual evaluation, did not converge or ended where the trial
# clamp departs from the identity. The surface Newton evaluates the condition
# at the clamped iterate and never warns: on a trial state the clamp is its
# purpose, on an accepted state it would be a failure.
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

# The mass leaving through the live boundary faces at the residual's last
# evaluation, per unit reference cross section, as `mass_loss_rate` sums it.
function _vbo_outflow(ds, c, topo, A0)
    out = 0.0
    for tag in Pyrolysis.Meshes.boundary_tags(topo),
        f0 in Pyrolysis.Meshes.boundary_faces(topo, tag)

        f = Pyrolysis.Meshes.live_face(topo, ds.active, f0)
        f == 0 || (out += sum(view(c.F_mass, :, f)) * ds.geometry.face_area[f])
    end
    return out / A0
end

# The mass ledger over the accepted steps: the outflow integrated
# over each step by three-point Gauss–Legendre on the step's interpolant,
# against the drop of the live cells' masses less the masses deactivation
# vented. A continuous callback whose condition never crosses zero reads each
# step before its events: the integral over the step and, when the threshold
# margin or the layer times' activation margin crosses inside it, over the
# step up to the first crossing, where a located event ends it.
# `commit!(integ)`, run after the step's events, adds the integral of the step
# as it ended; a step it cannot account for leaves the integral `NaN`.
mutable struct _VBOLedger
    outflow::Float64
    t_done::Float64
    on::Bool
    t::Float64
    tprev::Float64
    step::Float64
    crossing::Bool
    to_crossing::Float64
end

const _VBO_GAUSS = ((-sqrt(3 / 5), 5 / 9), (0.0, 8 / 9), (sqrt(3 / 5), 5 / 9))

function _vbo_ledger(cp)
    def, ws = cp.def, cp.ws
    topo, A0 = def.setup.topology, def.setup.reference.cross_section
    du = similar(cp.u0)
    rate(u, t) = begin
        Pyrolysis.Assembly.residual!(du, u, def, t, ws)
        _vbo_outflow(Pyrolysis.Assembly.buffers(ws, def, Float64)..., topo, A0)
    end
    integral(integ, a, b) = begin
        s = 0.0
        for (x, w) in _VBO_GAUSS
            τ = (a + b) / 2 + (b - a) / 2 * x
            s += w * rate(integ(τ), τ)
        end
        (b - a) / 2 * s
    end
    cells, layout = cp.deactivation, def.layout
    margins = Function[]
    cells === nothing ||
        push!(margins, u -> Pyrolysis.Deactivation.threshold_margin(cells, layout, u))
    corr = hasproperty(def.closure, :correction) ? def.closure.correction : nothing
    if corr isa Pyrolysis.Closures.LayerCorrection &&
       Pyrolysis.Closures.estimates_layer_times(corr)
        reference = def.setup.reference
        push!(
            margins,
            u -> Pyrolysis.Closures.layer_time_margin(corr, layout, reference, u),
        )
    end
    L = _VBOLedger(0.0, cp.problem.tspan[1], false, NaN, NaN, NaN, false, NaN)
    reader = ContinuousCallback(
        (u, t, integ) -> begin
            if t == integ.t > integ.tprev && !(L.on && integ.tprev == L.tprev)
                L.on, L.t, L.tprev = true, t, integ.tprev
                L.step = integral(integ, integ.tprev, t)
                cut = Inf
                for margin in margins
                    (margin(u) <= 0 && margin(integ.uprev) > 0) || continue
                    a, b = integ.tprev, t
                    for _ = 1:60
                        mid = (a + b) / 2
                        margin(integ(mid)) > 0 ? (a = mid) : (b = mid)
                    end
                    cut = min(cut, b)
                end
                L.crossing = isfinite(cut)
                L.crossing && (L.to_crossing = integral(integ, integ.tprev, cut))
            end
            return 1.0
        end,
        integ -> nothing;
        save_positions = (false, false),
        interp_points = 0,
    )
    commit!(integ) = begin
        if L.on && integ.tprev == L.t_done && integ.t == L.t
            L.outflow += L.step
        elseif L.on && integ.tprev == L.t_done && integ.t < L.t && L.crossing
            L.outflow += L.to_crossing
        else
            L.outflow = NaN
        end
        L.t_done, L.on = integ.t, false
    end
    return L, reader, commit!
end

# The solve with the step recorder, which evaluates the residual at every
# accepted state and reads the live cells' temperatures and the surface
# balances, and the mass ledger. `instrument(cp)` returns a further callback
# of the set-up problem for an instrumented measurement (`nothing` for none);
# `min_thickness` and, when given, `located_events` go to the integration,
# the rest to the setup.
function _vbo_recorded_solve(
    problem,
    rec::_VBORecord;
    min_thickness = 1e-6,
    located_events = nothing,
    instrument = cp -> nothing,
    kwargs...,
)
    cp = Pyrolysis.Integration.setup(
        problem;
        energy_state = VERIFICATION_ENERGY_STATE,
        kwargs...,
    )
    def, ws = cp.def, cp.ws
    topo = def.setup.topology
    thermal =
        Pyrolysis.Kernels.find_kernel(Pyrolysis.Kernels.ThermalBoundaryKernel, def.kernels)
    du = similar(cp.u0)
    live_mass(u, t) = begin
        Pyrolysis.Assembly.residual!(du, u, def, t, ws)
        sum(Pyrolysis.Assembly.buffers(ws, def, Float64)[1].m)
    end
    m_first = live_mass(cp.u0, cp.problem.tspan[1])
    m_last = m_first
    ledger, reader, commit! = _vbo_ledger(cp)
    observe!(integ) = begin
        commit!(integ)
        Pyrolysis.Assembly.residual!(du, integ.u, def, integ.t, ws)
        ds, _ = Pyrolysis.Assembly.buffers(ws, def, Float64)
        rec.T_max =
            max(rec.T_max, maximum(ds.T[i] for i in eachindex(ds.T) if ds.active[i]))
        rec.clamped += _vbo_clamped(thermal, ds, topo)
        m_last = sum(ds.m)
    end
    cp.deactivation === nothing || (rec.threshold = cp.deactivation.threshold)
    extra = instrument(cp)
    recorder = CallbackSet(reader, _vbo_recorder(observe!, rec))
    callback = extra === nothing ? recorder : CallbackSet(recorder, extra)
    located = located_events === nothing ? (;) : (; located_events)
    sol = Pyrolysis.Integration.integrate(cp; callback, min_thickness, located...)
    A0 = def.setup.reference.cross_section
    vented = sum((sum(e.lost_mass) for e in sol.events); init = 0.0)
    drop = (m_first - m_last) / A0
    rec.ledger = abs(drop - (ledger.outflow + vented / A0)) / drop
    return sol
end

# The case's solve. A moving mesh with deactivation gets the equidistribution
# correction, the layer interface Lagrangian, by `setup`'s installation rule.
function _vbo_solve_fn(id::String)
    return (problem; kwargs...) -> begin
        rec = _VBORecord()
        sol = _vbo_recorded_solve(problem, rec; kwargs...)
        VBO_STEP_STATS[id] = rec
        return sol
    end
end

# The peak and its time are the main peak's: over the saved times before the
# solve's last deactivation event, whose merge releases the layer's residue in
# a burst judged by refinement instead.
# `l_kw` is the substrate's thickness, the live thickness once the layer has
# burned out. The mass ledger is the case's recorder's (`id`; `NaN` without
# one).
function _vbo_measures(sol; l_kw = VBO_L_KW, id = nothing)
    mlr = sol.extras.gas_generation_rate .* 1000           # g m⁻² s⁻¹
    deactivations = sol.extras.deactivation_times
    before =
        isempty(deactivations) ? length(sol.t) :
        something(findlast(<(maximum(deactivations)), sol.t), 1)
    ipk = argmax(view(mlr, 1:before))
    th = sol.extras.thickness_history
    i_bo = findfirst(<(l_kw + 0.02e-3), th)
    t_bo = i_bo === nothing ? NaN : sol.t[i_bo]
    m = sol.extras.total_mass
    mass_loss = 100 * (1 - m[end] / m[1])
    ledger = id === nothing ? NaN : VBO_STEP_STATS[id].ledger
    return (;
        mlr,
        peak = mlr[ipk],
        t_peak = sol.t[ipk],
        t_bo,
        mass_loss,
        ledger,
        L_end = th[end],
        l_kw,
    )
end

# The mass loss of a deck whose layer of thickness `l_pmma` gasifies
# completely off the substrate of thickness `l_kw` [percent].
_vbo_mass_loss_exact(l_pmma, l_kw) =
    100 * VBO_RHO_PMMA * l_pmma / (VBO_RHO_PMMA * l_pmma + VBO_RHO_KW * l_kw)

# Maximum live-cell temperature at each output time (deactivated cells read
# NaN).
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

# The mass ledger's comment.
const VBO_LEDGER_COMMENT = "total-mass drop vs ∫ ṁ″ dt over the accepted steps (three-point Gauss–Legendre per step) and the mass deactivation vented"

# The exact and conservation QoIs shared by the burnout cases.
# The ledger closes to the integrator's tolerance: on the moving-mesh cases
# at their CI resolution at most 5.2e-7 in the temperature form and 3.2e-6 in
# the enthalpy form (V-BO-08's).
function _vbo_common_qois(
    m;
    id,
    mass_loss_exact = VBO_MASS_LOSS_EXACT,
    surface_warnings = true,
    ledger_tol = 3e-5,
    ledger_comment = VBO_LEDGER_COMMENT,
)
    out = QoIResult[]
    push!(
        out,
        QoIResult(
            name = "complete burnout: live thickness at t_end is the substrate",
            value = m.L_end,
            exact = m.l_kw,
            error = abs(m.L_end - m.l_kw),
            tol = 1e-5,
            comment = "every PMMA cell deactivated; $(round(1e3 * m.l_kw; digits = 2)) ± 0.01 mm",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "mass loss of the gasified PMMA layer",
            value = m.mass_loss,
            exact = mass_loss_exact,
            error = abs(m.mass_loss - mass_loss_exact),
            tol = 0.3,
            comment = "percent of the initial mass; exact = the layer's initial mass / m₀",
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

# The surface energy balance's QoI: the Newton clamps its trial iterates
# smoothly and never warns, so no accepted state may hold a balance that did
# not converge or sits on the clamp.
function _vbo_surface_qoi(id)
    n = get(() -> _VBORecord(-1.0, -1.0, NaN, -1, NaN, NaN), VBO_STEP_STATS, id).clamped
    return QoIResult(
        name = "clamped surface balances on accepted states",
        value = Float64(n),
        exact = 0.0,
        error = Float64(abs(n)),
        tol = 0.5,
        comment = "boundary temperature solve converged inside the property bracket at every accepted step",
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

# The deactivation invocations of the PMMA layer, with the cells they took
# and the largest merged fraction, reported and not asserted: the burnout
# path's criterion is its accuracy under refinement (test/gates/burnout_envelope.jl),
# and one invocation is its target.
function _vbo_invocation_qoi(sol, id)
    ε_m = VBO_STEP_STATS[id].threshold
    inv = _vbo_invocations(sol, ε_m)
    largest = maximum(v.largest for v in inv; init = 0.0)
    return QoIResult(
        name = "deactivation invocations for the PMMA layer (reported)",
        value = Float64(length(inv)),
        exact = 1.0,
        error = abs(length(inv) - 1.0),
        tol = Inf,
        comment = "$(sum(v.cells for v in inv; init = 0)) cells; largest merged fraction $(round(largest; sigdigits = 3)) ($(round(largest / ε_m; sigdigits = 3)) ε_m); reported, not asserted",
    )
end

# The report figures of a burnout case: its gas generation rate and live
# thickness through the run.
_vbo_history_curves(id, sol) = [
    _vbo_history_curve(
        id,
        :mlr_history,
        "PMMA gasification rate through burnout",
        "MLR (g m⁻² s⁻¹)",
        "gas generation rate",
        sol.t,
        sol.extras.gas_generation_rate .* 1000,
    ),
    _vbo_history_curve(
        id,
        :thickness_history,
        "live thickness",
        "thickness (m)",
        "live thickness",
        sol.t,
        sol.extras.thickness_history,
    ),
]

# V-BO-01's solve keywords, shared by the burnout cases that vary its deck.
const VBO_SOLVE_KWARGS = (
    use_ale = true,
    deactivation = true,
    min_thickness = 1e-3,
    radiation_model = BEER_LAMBERT,
    abstol = 1e-11,
    reltol = 1e-9,
    pressure_tolerance = 1e-4,
)

function _verification_case_vbo01()
    id = "V-BO-01"
    function qois(sol)
        m = _vbo_measures(sol; id)
        out = _vbo_common_qois(m; id)
        push!(out, _vbo_invocation_qoi(sol, id))
        return out
    end

    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout under ALE (layer deactivation)",
        tags = [:ale, :burnout, :depletion, :beer_lambert, :coupled, :performance],
        references = [
            "MaCFP-2 PMMA gasification case (NIST), 50 kW/m²",
            "Gpyro V2 PMMA deck property fits",
        ],
        build = n -> _vbo_problem(n),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = VBO_SOLVE_KWARGS,
        qois = qois,
        # The tight pressure/ledger controls cost about 111 s for a warmed
        # enthalpy run at 60 cells (8690 accepted steps, 1154 rejections).
        ci = CIConfig(n_cells = 60, budget_s = 120.0),
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
