# =============================================================================
# V-BO-03 — PMMA burnout on a transparent substrate (transmitted-flux guard)
# =============================================================================
#
# The V-BO-01 deck under ALE with the Kaowool absorption coefficient set to
# zero, so the substrate neither absorbs nor emits. Radiation that passes the
# PMMA layer, the heater's and the PMMA's own, crosses the substrate and
# leaves through the back face without heating it, so the share of the
# radiation entering the exposed face that leaves through the back rises
# toward one as the PMMA thins. The mesh moves by the equidistribution
# correction as in V-BO-01. The run must complete with the surface
# temperature converged at every accepted step, and the radiation leaving
# through the back face (solution.extras.radiation_back) must reach more than
# 1 % of the radiation entering (solution.extras.radiation_absorbed).

# The end the burnout measurements run V-BO-03's deck to, past its layer's
# last invocation (898, 817 and 856 s at 30, 60 and 120 cells).
const VBO03_T_EXTENDED = 1000.0

# The share of the radiation entering the exposed face that leaves through the
# back face, at every saved time.
_vbo03_back_share(sol) = sol.extras.radiation_back ./ sol.extras.radiation_absorbed

function _verification_case_vbo03()
    id = "V-BO-03"
    function qois(sol)
        m = _vbo_measures(sol)
        out = QoIResult[_vbo_surface_qoi(id)]
        leaving = maximum(_vbo03_back_share(sol))
        push!(
            out,
            QoIResult(
                name = "radiation leaving through the back reported (series peak above 1 % of the radiation entering)",
                value = leaving,
                exact = 1.0,
                error = leaving > 0.01 ? 0.0 : 1.0,
                tol = 0.5,
                comment = "largest share of the radiation entering the exposed face that leaves through the transparent back face",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "peak share of the radiation entering that leaves through the back",
                value = leaving,
                exact = 1.0,
                error = 1.0 - leaving,
                tol = 0.3,
                comment = "asymptote 1 once the PMMA is gone, since the substrate neither absorbs nor emits; at t_end the PMMA cells next to the substrate still hold about 1 % of their mass and absorb 8 % of it, and 3 % leaves through the exposed face (observed 0.889)",
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
                comment = "percent; the radiation crossing the substrate does not heat it",
            ),
        )
        return out
    end
    return VerificationCase(
        id = id,
        title = "PMMA burnout on a transparent substrate (transmitted-flux guard)",
        tags = [:ale, :burnout, :depletion, :radiative_transfer],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m²"],
        build = n -> _vbo_problem(n; substrate_alpha = 0.0),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = (
            use_ale = true,
            deactivation = true,
            min_thickness = 1e-3,
            radiation_model = RADIATIVE_TRANSFER,
            abstol = 1e-8,
            reltol = 1e-6,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 60, budget_s = 10.0),
        report = (
            n_cells_ladder = [30, 60, 120],
            curves = sol -> [
                _vbo_history_curve(
                    id,
                    :back_history,
                    "radiation leaving through the back face",
                    "share of the radiation entering the exposed face",
                    "leaving through the back",
                    sol.t,
                    _vbo03_back_share(sol),
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
