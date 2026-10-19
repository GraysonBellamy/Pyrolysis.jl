# =============================================================================
# V-BO-02 — PMMA-on-Kaowool burnout on a fixed mesh under surface absorption
# =============================================================================
#
# The V-BO-01 deck on a fixed mesh, with the incident flux absorbed at the
# surface. A fixed mesh keeps every cell's volume, so a gasifying cell empties
# instead of thinning, and only the top live cell parks: once its condensed
# concentration falls below 5 % of its initial value it merges into the cell
# beneath, which becomes the surface. With the flux absorbed at the surface,
# the surface cell is the hottest and empties first, so the layer burns top
# down and no cell below the surface empties before it.
#
# The run must park all 60 PMMA cells, lose exactly the mass of the PMMA
# layer, end on the bare substrate, close the mass ledger, keep the maximum
# live-cell temperature below 1000 K, and take no step below the post-park
# reset. The peak mass loss rate is informational: surface absorption is a
# different problem from the in-depth reference of V-BO-01.

# Peak of the CI run (0.098 mm cells).
const VBO02_PEAK_REF = 24.62    # g m⁻² s⁻¹, at 429.5 s

function _verification_case_vbo02()
    id = "V-BO-02"
    function qois(sol)
        m = _vbo_measures(sol)
        out = _vbo_common_qois(
            m;
            id,
            peak_ref = VBO02_PEAK_REF,
            peak_tol = 0.25,
            peak_comment = "g m⁻² s⁻¹; surface absorption on a fixed mesh, a different problem from the in-depth reference (informational)",
        )
        st = sol.raw_solution.stats
        rec = VBO_STEP_STATS[id]
        n_park = length(sol.extras.park_times)
        n_pmma = count(>(0.0), view(sol.ξ, 1, :, 1))
        push!(
            out,
            QoIResult(
                name = "park events",
                value = Float64(n_park),
                exact = Float64(n_pmma),
                error = abs(n_park - n_pmma),
                tol = 0.5,
                comment = "every PMMA cell parks as the top live cell, top down",
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
                exact = Pyrolysis.Solver.DT_AFTER_PARK,
                error = max(0.0, 1.0 - rec.dt_min / Pyrolysis.Solver.DT_AFTER_PARK),
                tol = 1e-3,
                comment = "s; no step below the post-park reset (observed 1e-6 s)",
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
                comment = "budget at the CI resolution (60 PMMA cells); observed 56, about one per park interval",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "accepted steps",
                value = Float64(st.naccept),
                exact = 0.0,
                error = Float64(st.naccept),
                tol = VERIFICATION_CORE == "old" ? 2000.0 :
                      VERIFICATION_ENERGY_STATE === :enthalpy ? 4000.0 : 3000.0,
                comment = "budget at the CI resolution (60 PMMA cells); observed 1572 on the old core, about 26 per park interval, and on the new one 2354 (temperature form) and 3163 (enthalpy form)",
            ),
        )
        return out
    end
    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout on a fixed mesh under surface absorption (top-cell parking)",
        tags = [:eulerian, :burnout, :depletion, :surface_absorption, :coupled],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m²"],
        build = n -> _vbo_problem(n; substrate_alpha = 0.0, pmma_alpha = 0.0),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = (
            use_ale = false,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
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
