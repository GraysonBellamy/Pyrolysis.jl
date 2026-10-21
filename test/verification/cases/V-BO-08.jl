# =============================================================================
# V-BO-08 — localized burnout of a thin layer on a deep substrate
# =============================================================================
#
# The V-BO-01 deck with a PMMA layer of ten cells at V-BO-01's pitch (about
# 1 mm) on the 12 mm Kaowool. The layer burns out at about 70 s and the
# substrate heats on beneath the bare surface to 150 s. The run must lose
# exactly the layer's mass and end on the bare Kaowool. On the finite-volume
# core the layer
# deactivates in one invocation and the correction stays where the layer is:
# the Kaowool, which cannot vanish, keeps its mesh with the solid, so the
# correction velocity is exactly zero at every Kaowool node at every saved
# time.

const VBO08_N_PMMA = 10
const VBO08_L_PMMA = VBO08_N_PMMA * VBO_L_PMMA / 60
const VBO08_T_END = 150.0

function _verification_case_vbo08()
    id = "V-BO-08"
    function qois(sol)
        l_kw = first(sol.extras.thickness_history) - VBO08_L_PMMA
        m = _vbo_measures(sol; l_kw, id)
        out = _vbo_common_qois(
            m;
            id,
            mass_loss_exact = _vbo_mass_loss_exact(VBO08_L_PMMA, l_kw),
        )
        VERIFICATION_CORE == "new" || return out
        push!(out, _vbo_invocation_qoi(sol, id))
        # The substrate's nodes: those at or below its initial top.
        n_kw = count(<=(l_kw * (1 + 1e-9)), view(sol.z, :, 1))
        c = view(sol.correction_velocity, 1:n_kw, :)
        push!(
            out,
            QoIResult(
                name = "the correction is zero at every Kaowool node",
                value = maximum(abs, c),
                exact = 0.0,
                error = maximum(abs, c),
                tol = eps(),
                comment = "m/s over the $(n_kw) substrate nodes and every saved time",
            ),
        )
        return out
    end
    return VerificationCase(
        id = id,
        title = "Localized burnout of a thin PMMA layer on a deep substrate",
        tags = [:ale, :burnout, :depletion, :beer_lambert, :coupled],
        references = ["MaCFP-2 PMMA gasification case (NIST), 50 kW/m², a 1 mm layer"],
        build = n -> _vbo_problem(n; l_pmma = VBO08_L_PMMA, t_end = VBO08_T_END),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = VBO_SOLVE_KWARGS,
        qois = qois,
        ci = CIConfig(n_cells = VBO08_N_PMMA, budget_s = 20.0),
        report = (curves = sol -> _vbo_history_curves(id, sol),),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo08)
