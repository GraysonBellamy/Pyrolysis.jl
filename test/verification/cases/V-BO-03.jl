# =============================================================================
# V-BO-03 — PMMA burnout on a transparent substrate (transmitted-flux guard)
# =============================================================================
#
# The V-BO-01 deck under ALE with the Kaowool absorption coefficient set to
# zero. Radiation that passes the PMMA layer crosses the substrate and leaves
# through the back face without heating it, so the share of the incident flux
# that leaves the slab rises toward one as the PMMA thins. The run must
# complete, report the flux that leaves through the back face
# (solution.extras.radiation_transmitted) and warn about it once, with the
# surface temperature converged at every evaluation.

function _verification_case_vbo03()
    id = "V-BO-03"
    function qois(sol)
        m = _vbo_measures(sol)
        out = QoIResult[]
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
        n_tr = get(VBO_TRANSMITTED_WARNINGS, id, -1)
        push!(
            out,
            QoIResult(
                name = "transmitted-flux warning fired",
                value = Float64(n_tr),
                exact = 1.0,
                error = abs(n_tr - 1.0),
                tol = 0.5,
                comment = "radiation leaving through the transparent back face is reported",
            ),
        )
        tr = sol.extras.radiation_transmitted
        push!(
            out,
            QoIResult(
                name = "peak transmitted fraction of the incident flux",
                value = maximum(tr) / (0.96 * 50e3),
                exact = 1.0,
                error = 1.0 - maximum(tr) / (0.96 * 50e3),
                tol = 0.3,
                comment = "asymptote 1 once the PMMA is gone; the transmitted flux does not heat the substrate, so a residual PMMA layer remains at t_end and absorbs about 20 %",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "mass loss (informational)",
                value = m.mass_loss,
                exact = VBO_MASS_LOSS_EXACT,
                error = abs(m.mass_loss - VBO_MASS_LOSS_EXACT),
                tol = 10.0,
                comment = "percent; the transmitted flux does not heat the substrate",
            ),
        )
        return out
    end
    return VerificationCase(
        id = id,
        title = "PMMA burnout on a transparent substrate (transmitted-flux guard)",
        tags = [:ale, :burnout, :depletion, :beer_lambert],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m²"],
        build = n -> _vbo_problem(n; substrate_alpha = 0.0),
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
            curves = sol -> [
                _vbo_history_curve(
                    id,
                    :transmitted_history,
                    "flux leaving through the back face",
                    "transmitted / absorbed incident flux",
                    "transmitted fraction",
                    sol.t,
                    sol.extras.radiation_transmitted ./ (0.96 * 50e3),
                ),
                _vbo_history_curve(
                    id,
                    :mlr_history,
                    "PMMA gasification rate on a transparent substrate",
                    "MLR (g m⁻² s⁻¹)",
                    "gas generation rate",
                    sol.t,
                    sol.extras.gas_generation_rate .* 1000,
                ),
            ],
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo03)
