# =============================================================================
# V-BO-02 — PMMA-on-Kaowool burnout on a fixed mesh under surface absorption
# =============================================================================
#
# The V-BO-01 deck on a fixed mesh, with the incident flux absorbed at the
# surface. A fixed mesh keeps every cell's volume, so a gasifying cell empties
# instead of thinning, and only the top live cell deactivates: once its
# condensed mass falls below ε_m (1 %) of its initial value it merges into the
# cell beneath, which becomes the surface. With the flux absorbed at the
# surface, the surface cell is the hottest and empties first, so the layer
# burns top down and no cell below the surface empties before it.
#
# The run must deactivate all 60 PMMA cells, lose exactly the mass of the
# PMMA layer, end on the bare substrate, close the mass ledger, keep the
# maximum live-cell temperature below 1000 K, and take no step below the
# restart step after an event. The peak mass loss rate is informational:
# surface absorption is a different problem from the in-depth reference of
# V-BO-01.

# Peak of the CI run (0.098 mm cells).
const VBO02_PEAK_REF = 23.76    # g m⁻² s⁻¹, at 438.5 s

function _verification_case_vbo02()
    id = "V-BO-02"
    function qois(sol)
        m = _vbo_measures(sol; id)
        out = _vbo_common_qois(
            m;
            id,
            peak_ref = VBO02_PEAK_REF,
            peak_tol = 0.25,
            peak_comment = "g m⁻² s⁻¹; surface absorption on a fixed mesh, a different problem from the in-depth reference (informational)",
            ledger_tol = 3e-4,
            ledger_comment = VBO_LEDGER_COMMENT *
                             "; each deactivation leaves a gas transient on the fixed mesh, so this is the time error at reltol 1e-6 (observed 6.1e-5 and 6.2e-5 in the two energy forms at 60 cells, 1.2e-4 at 120, 1.9e-8 at reltol 1e-8)",
        )
        st = sol.raw_solution.stats
        rec = VBO_STEP_STATS[id]
        n_events = length(sol.extras.deactivation_times)
        n_pmma = count(>(0.0), view(sol.ξ, 1, :, 1))
        push!(
            out,
            QoIResult(
                name = "deactivation events",
                value = Float64(n_events),
                exact = Float64(n_pmma),
                error = abs(n_events - n_pmma),
                tol = 0.5,
                comment = "every PMMA cell deactivates as the top live cell, top down",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "maximum live-cell temperature",
                value = rec.T_max,
                error = rec.T_max,
                tol = 1000.0,
                comment = "K over every accepted step (observed 915); the surface cell is the hottest and empties first, so no interior void forms",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "smallest accepted step",
                value = rec.dt_min,
                exact = Pyrolysis.Integration.DT_AFTER_EVENT,
                error = max(0.0, 1.0 - rec.dt_min / Pyrolysis.Integration.DT_AFTER_EVENT),
                tol = 1e-3,
                comment = "s; no step below the restart step after an event, DT_AFTER_EVENT (observed 5.0e-5 s in both energy forms)",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "Newton failures",
                value = Float64(st.nnonlinconvfail),
                exact = 0.0,
                error = Float64(st.nnonlinconvfail),
                tol = 100.0,
                comment = "budget at the CI resolution (60 PMMA cells); observed 0: the default integrator, Rodas5P, is a Rosenbrock method and solves no nonlinear system",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "accepted steps",
                value = Float64(st.naccept),
                exact = 0.0,
                error = Float64(st.naccept),
                tol = VERIFICATION_ENERGY_STATE === :enthalpy ? 4000.0 : 3000.0,
                comment = "budget at the CI resolution (60 PMMA cells); observed 1978 (temperature form) and 2797 (enthalpy form)",
            ),
        )
        return out
    end
    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout on a fixed mesh under surface absorption (top-cell deactivation)",
        tags = [:eulerian, :burnout, :depletion, :surface_absorption, :coupled],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m²"],
        build = n -> _vbo_problem(n; substrate_alpha = 0.0, pmma_alpha = 0.0),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = (
            use_ale = false,
            deactivation = true,
            radiation_model = SURFACE_ABSORPTION,
            abstol = 1e-8,
            reltol = 1e-6,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 60, budget_s = 20.0),
        report = (
            n_cells_ladder = [30, 60, 120],
            curves = sol -> [
                _vbo_history_curve(
                    id,
                    :mlr_history,
                    "PMMA gasification rate, top-down burnout",
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
                _vbo_history_curve(
                    id,
                    :max_temperature_history,
                    "maximum live-cell temperature",
                    "T (K)",
                    "maximum live-cell temperature",
                    sol.t,
                    _vbo_max_live_T(sol),
                ),
            ],
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo02)
