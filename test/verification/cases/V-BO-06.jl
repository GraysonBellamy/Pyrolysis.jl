# =============================================================================
# V-BO-06 — V-BO-01 on a mesh graded toward the exposed surface
# =============================================================================
#
# The V-BO-01 deck with its PMMA cells graded geometrically toward the exposed
# surface, each cell 1/1.02 of the one beneath it at 60 cells, so the surface
# cell is a third of the cell at the PMMA–Kaowool interface; at other
# resolutions the ratio per cell keeps that third, so the cell ladder refines
# one mesh. The Kaowool takes the interface cell's pitch. The burnout
# quantities do not depend on the mesh, so V-BO-01's references and
# tolerances hold. The equidistribution correction starts from cells of
# unequal thickness, and the invocations that deactivate the layer are
# reported.

const VBO06_GRADING = 1.02
const VBO06_CELLS = 60

# The ratio per cell at `n` PMMA cells: the surface cell keeps the 60-cell
# mesh's share of the interface cell.
_vbo06_grading(n) = VBO06_GRADING^((VBO06_CELLS - 1) / (n - 1))

function _verification_case_vbo06()
    id = "V-BO-06"
    n_default = VBO06_CELLS
    function qois(sol)
        l_kw = first(sol.extras.thickness_history) - VBO_L_PMMA
        m = _vbo_measures(sol; l_kw, id)
        out = _vbo_common_qois(
            m;
            id,
            peak_tol = 0.03,
            mass_loss_exact = _vbo_mass_loss_exact(VBO_L_PMMA, l_kw),
        )
        append!(out, _vbo_timing_qois(m))
        push!(out, _vbo_invocation_qoi(sol, id))
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
        build = n -> _vbo_problem(n; grading = _vbo06_grading(n)),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = VBO_SOLVE_KWARGS,
        qois = qois,
        ci = CIConfig(n_cells = n_default, budget_s = 20.0),
        report = (curves = sol -> _vbo_history_curves(id, sol),),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo06)
