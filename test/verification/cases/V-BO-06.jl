# =============================================================================
# V-BO-06 — V-BO-01 on a mesh graded toward the exposed surface
# =============================================================================
#
# The V-BO-01 deck with its PMMA cells graded geometrically toward the exposed
# surface, each cell 1/1.02 of the one beneath it, so the surface cell is a
# third of the cell at the PMMA–Kaowool interface; the Kaowool takes the
# interface cell's pitch. The burnout quantities do not depend on the mesh, so
# V-BO-01's references and tolerances hold. On the finite-volume core the
# equidistribution correction starts from cells of unequal thickness, and the
# layer must still deactivate in one invocation. The old core's ledger is
# bounded at 5e-3: its first parks, of the thin surface cells, overstate the
# outflow by 0.4 % of the mass lost whatever the output interval (0.08 % on
# V-BO-01's uniform mesh), where the finite-volume core closes it to 1e-4.

const VBO06_GRADING = 1.02

function _verification_case_vbo06()
    id = "V-BO-06"
    n_default = 60
    function qois(sol)
        l_kw = first(sol.extras.thickness_history) - VBO_L_PMMA
        m = _vbo_measures(sol; l_kw)
        old = VERIFICATION_CORE == "old"
        out = _vbo_common_qois(
            m;
            id,
            peak_tol = 0.03,
            mass_loss_exact = _vbo_mass_loss_exact(VBO_L_PMMA, l_kw),
            ledger_tol = old ? 5e-3 : 1e-3,
            ledger_comment = "total-mass drop vs ∫ ṁ″ dt over the run" * (
                old ?
                "; the old core's first parks overstate the outflow (observed 3.8e-3)" :
                ""
            ),
        )
        append!(out, _vbo_timing_qois(m))
        VERIFICATION_CORE == "new" && push!(out, _vbo_invocation_qoi(sol, id))
        return out
    end
    return VerificationCase(
        id = id,
        title = "PMMA-on-Kaowool burnout on a mesh graded toward the exposed surface",
        tags = [:ale, :burnout, :depletion, :beer_lambert, :coupled],
        references = [
            "MaCFP-2 PMMA gasification case (NIST), 50 kW/m²",
            "Gpyro V2 PMMA deck property fits",
        ],
        build = n -> _vbo_problem(n; grading = VBO06_GRADING),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = VBO_SOLVE_KWARGS,
        qois = qois,
        ci = CIConfig(n_cells = n_default, budget_s = 20.0),
        report = (curves = sol -> _vbo_history_curves(id, sol),),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo06)
