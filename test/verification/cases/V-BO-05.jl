# =============================================================================
# V-BO-05 — fixed-mesh burnout under radiative transfer: the layer gasifies from within
# =============================================================================
#
# The V-BO-01 deck, radiative transfer included, on a fixed mesh. The cells
# absorb and emit through the depth while the exposed face also loses heat by
# convection, so the temperature maximum of the PMMA lies below the surface
# and the PMMA beneath the surface empties nearly uniformly, the exposed cell
# slowest.
#
# A cell deactivates when its reactive inventory falls below ε_m of its
# reference. A fixed mesh keeps every cell's volume, so a cell emptied beneath
# one that holds its material is a cavity and would end the solve with outcome
# `:interior`; but the cells above the first crossing here hold only 1.0 to 1.4
# thresholds, below the band, so the crossing cell leaves together with them,
# top down into the cell beneath it, and the layer burns through from within.
# At the CI resolution, 60 PMMA cells, the cell beneath the exposed one crosses
# first at 157.5 s with the exposed cell holding 1.2 thresholds, the two go in
# one invocation, and the layer is gone by 372 s in 23 invocations: every PMMA
# cell deactivates, the run loses exactly the layer's mass, ends on the bare
# substrate and closes the mass ledger, with the maximum live-cell temperature
# below 1000 K. The histories and the profiles at the last output before the
# first event, the skin at the top over the emptying layer, are reported.

function _vbo05_report_curves(sol)
    id = "V-BO-05"
    t_first = isempty(sol.events) ? sol.t[end] : first(sol.events).t
    k = something(findlast(<(t_first), sol.t), length(sol.t))
    z_mm = sol.cell_centers[:, k] .* 1e3
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
            "PMMA gasification rate",
            "MLR (g m⁻² s⁻¹)",
            "gas generation rate",
            sol.t,
            sol.extras.gas_generation_rate .* 1000,
        ),
        _vbo_history_curve(
            id,
            :max_temperature_history,
            "maximum live-cell temperature",
            "T (K)",
            "maximum live-cell temperature",
            sol.t,
            _vbo_max_live_T(sol),
        ),
        profile(
            :pmma_profile,
            "PMMA concentration before the first event",
            "ξ_PMMA (kg m⁻³)",
            "PMMA concentration",
            view(sol.ξ, 1, :, k),
        ),
        profile(
            :temperature_profile,
            "temperature before the first event",
            "T (K)",
            "temperature",
            view(sol.T, :, k),
        ),
    ]
end

# The QoIs: complete burnout with the layer's mass and the ledger (shared with
# V-BO-01), every PMMA cell deactivated in fewer invocations than cells (the
# runs take the spent cells above a crossing), no cavity, and the temperature
# and work bounds.
function _vbo05_qois(sol)
    id = "V-BO-05"
    m = _vbo_measures(sol; id)
    out = _vbo_common_qois(
        m;
        id,
        ledger_tol = 3e-4,
        ledger_comment = VBO_LEDGER_COMMENT *
                         "; the fixed mesh vents each dead cell's gas to its event",
    )
    st = sol.raw_solution.stats
    rec = VBO_STEP_STATS[id]
    n_events = length(sol.extras.deactivation_times)
    n_pmma = count(>(0.0), view(sol.ξ, 1, :, 1))
    invocations = sol.statistics.invocations
    push!(
        out,
        QoIResult(
            name = "deactivation events",
            value = Float64(n_events),
            exact = Float64(n_pmma),
            error = abs(n_events - n_pmma),
            tol = 0.5,
            comment = "every PMMA cell deactivates; the layer gasifies from within and leaves with its spent skins",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "invocations fewer than events",
            value = Float64(invocations),
            error = invocations < n_events ? 0.0 : 1.0,
            tol = 0.5,
            comment = "at least one run took a spent cell above its crossing (23 invocations for 60 cells at the CI resolution; the first at 157.5 s)",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "outcome",
            value = sol.outcome === :interior ? 1.0 : 0.0,
            exact = 0.0,
            error = sol.outcome === :interior ? 1.0 : 0.0,
            tol = 0.5,
            comment = "no cavity: the solve ended with outcome :$(sol.outcome)",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "maximum live-cell temperature",
            value = rec.T_max,
            error = rec.T_max,
            tol = 1000.0,
            comment = "K over every accepted step (observed 927 at 60 cells)",
        ),
    )
    push!(
        out,
        QoIResult(
            name = "accepted steps",
            value = Float64(st.naccept),
            exact = 0.0,
            error = Float64(st.naccept),
            tol = (VERIFICATION_ENERGY_STATE === :enthalpy ? 12000.0 : 10000.0) *
                  max(1, n_pmma / 60),
            comment = "work budget, scaled by cell count; 1856 temperature steps at 60 cells",
        ),
    )
    return out
end

function _verification_case_vbo05()
    id = "V-BO-05"
    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout on a fixed mesh under radiative transfer: the layer gasifies from within",
        tags = [:eulerian, :burnout, :depletion, :radiative_transfer, :coupled],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m²"],
        build = n -> _vbo_problem(n),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = (
            use_ale = false,
            deactivation = true,
            radiation_model = RADIATIVE_TRANSFER,
            abstol = 1e-8,
            reltol = 1e-6,
        ),
        qois = _vbo05_qois,
        ci = CIConfig(n_cells = 60, budget_s = 20.0),
        report = (n_cells_ladder = [30, 60, 120], curves = _vbo05_report_curves),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo05)
