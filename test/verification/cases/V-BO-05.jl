# =============================================================================
# V-BO-05 — fixed-mesh burnout under in-depth absorption stops at the first interior void
# =============================================================================
#
# The V-BO-01 deck, in-depth absorption included, on a fixed mesh. The
# temperature maximum of the PMMA lies below the surface, so the surface cell
# reacts slowest: it persists as a low-density skin while the cells beneath it
# empty, and the layer gasifies from within. A void cell conducts with the gas
# conductivity and absorbs nothing, and nothing in the model carries heat
# across a void, so the flux transmitted through the void column would heat
# the substrate beneath it far outside the property bracket.
#
# The depletion callback stops a fixed-mesh solve when a live cell below the
# top live cell becomes a void (matrix volume fraction below 1e-4). The run
# must stop between 330 and 345 s with the interior-void warning, before any
# park, with the maximum live-cell temperature below 1000 K. The solve ends
# with retcode Terminated.
#
# On the finite-volume core a cell deactivates when its condensed mass falls
# below ε_m of its initial value, and on a fixed mesh a crossing cell beneath
# a live cell ends the solve `:interior`, without a warning. That criterion
# stops the run earlier, at 259.6 s; the stop time's tolerance is three times
# the CI resolution's distance from the 240-cell time (260.35, 259.79, 259.62
# and 259.57 s at 30, 60, 120 and 240 cells, both energy forms within 5 ms).

# Histories up to the stop, and the PMMA and temperature profiles at the last
# output before it: the skin at the top, the emptying layer beneath.
function _vbo05_report_curves(sol)
    id = "V-BO-05"
    k = length(sol.t)
    z_mm = cell_positions(sol) .* 1e3
    t_k = round(sol.t[k], digits = 1)
    profile(name, title, ylabel, label, y) = ReportCurve(
        name = name,
        title = "$(id): $(title) at t = $(t_k) s",
        xlabel = "height above the back face (mm)",
        ylabel = ylabel,
        series = [CurveSeries(label = label, x = z_mm, y = collect(Float64, y))],
    )
    return [
        _vbo_history_curve(
            id,
            :mlr_history,
            "PMMA gasification rate up to the stop",
            "MLR (g m⁻² s⁻¹)",
            "gas generation rate",
            sol.t,
            sol.extras.gas_generation_rate .* 1000,
        ),
        _vbo_history_curve(
            id,
            :max_temperature_history,
            "maximum live-cell temperature up to the stop",
            "T (K)",
            "maximum live-cell temperature",
            sol.t,
            _vbo_max_live_T(sol),
        ),
        profile(
            :pmma_profile,
            "PMMA concentration",
            "ξ_PMMA (kg m⁻³)",
            "PMMA concentration",
            view(sol.ξ, 1, :, k),
        ),
        profile(
            :temperature_profile,
            "temperature",
            "T (K)",
            "temperature",
            view(sol.T, :, k),
        ),
    ]
end

# The finite-volume core's QoIs: the stop at the first ε_m crossing beneath a
# live cell, its outcome, no deactivation before it, and the temperature bound.
function _vbo05_core_qois(sol, rec)
    t_stop = rec.t_last
    n_events = length(sol.events)
    return [
        QoIResult(
            name = "stop time at the first interior crossing of ε_m",
            value = t_stop,
            exact = 259.57,
            error = abs(t_stop - 259.57),
            tol = 0.7,
            comment = "s; the 240-cell time (259.79 at 60 cells): a cell under the skin falls below ε_m of its initial condensed mass",
        ),
        QoIResult(
            name = "outcome :interior",
            value = sol.outcome == :interior ? 1.0 : 0.0,
            exact = 1.0,
            error = sol.outcome == :interior ? 0.0 : 1.0,
            tol = 0.5,
            comment = "the solve ended with outcome :$(sol.outcome)",
        ),
        QoIResult(
            name = "deactivation events before the stop",
            value = Float64(n_events),
            exact = 0.0,
            error = Float64(n_events),
            tol = 0.5,
            comment = "the surface cell persists as a skin while the layer gasifies from within",
        ),
        QoIResult(
            name = "maximum live-cell temperature up to the stop",
            value = rec.T_max,
            error = rec.T_max,
            tol = 1000.0,
            comment = "K over every accepted step (observed 693)",
        ),
    ]
end

function _verification_case_vbo05()
    id = "V-BO-05"
    function qois(sol)
        rec = VBO_STEP_STATS[id]
        VERIFICATION_CORE == "new" && return _vbo05_core_qois(sol, rec)
        t_stop = rec.t_last
        n_void = get(VBO_VOID_WARNINGS, id, -1)
        n_park = length(sol.extras.park_times)
        return [
            QoIResult(
                name = "stop time at the first interior void",
                value = t_stop,
                exact = 337.5,
                error = abs(t_stop - 337.5),
                tol = 7.5,
                comment = "s; window 330–345 s (observed 337.0): the lowest PMMA cell reaches matrix volume fraction < MIXING_BLEND_LO under the skin",
            ),
            QoIResult(
                name = "interior-void warning",
                value = Float64(n_void),
                exact = 1.0,
                error = abs(n_void - 1.0),
                tol = 0.5,
                comment = "the stop names the cause and the ALE mesh",
            ),
            QoIResult(
                name = "park events before the stop",
                value = Float64(n_park),
                exact = 0.0,
                error = Float64(n_park),
                tol = 0.5,
                comment = "the surface cell persists as a skin while the layer gasifies from within",
            ),
            QoIResult(
                name = "maximum live-cell temperature up to the stop",
                value = rec.T_max,
                error = rec.T_max,
                tol = 1000.0,
                comment = "K over every accepted step (observed 912); the solve stops before the void column runs away",
            ),
        ]
    end
    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout on a fixed mesh under in-depth absorption stops at the first interior void",
        tags = [:eulerian, :burnout, :depletion, :beer_lambert, :coupled],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m²"],
        build = n -> _vbo_problem(n),
        solve_fn = _vbo_solve_fn(id),
        expected_retcode = :Terminated,
        solve_kwargs = (
            use_ale = false,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
            radiation_model = BEER_LAMBERT,
            abstol = 1e-8,
            reltol = 1e-6,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 60, budget_s = 10.0),
        report = (n_cells_ladder = [30, 60, 120], curves = _vbo05_report_curves),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo05)
