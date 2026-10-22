# =============================================================================
# V-BO-07 — burnout of a heterogeneous layer of mixed composition
# =============================================================================
#
# The V-BO-01 deck with a second solid mixed into the PMMA layer: the same
# properties at half the density (605 kg/m³) and a tenth of the
# pre-exponential factor, its mass fraction rising linearly from zero at the
# PMMA–Kaowool interface to one half at the exposed surface, each cell's
# concentrations filling its volume. The layer's density therefore varies
# across one activated region, which is where the equidistribution monitor's
# mass term acts. Both solids gasify completely, so the run loses exactly the
# layer's initial mass and ends on the bare Kaowool. The references of the
# peak, its time and the burnout time are the old core's at 240 PMMA cells
# (the peak 30.41, 30.54, 30.61 and 30.64 g m⁻² s⁻¹ at 30, 60, 120 and 240).
#
# On the finite-volume core the layer deactivates in one invocation. The burn
# runs down from the surface: the upper cells lose their PMMA first and keep
# the slow solid, and the PMMA beneath burns last and fast. When the first cell
# crosses ε_m the others still hold up to about 6 % of their initial mass,
# above the 5 ε_m V-BO-01 bounds, with or without the monitor's mass term, so
# the case asserts the one invocation and reports that fraction.

const VBO07_RHO_LIGHT = VBO_RHO_PMMA / 2
const VBO07_Y_SURFACE = 0.5
const VBO07_PEAK_REF = 30.64     # g m⁻² s⁻¹
const VBO07_T_PEAK_REF = 288.0   # s
const VBO07_T_BO_REF = 325.0     # s

# The light solid's mass fraction at height `z` in a layer from `z_if` up.
_vbo07_fraction(z, z_if) = VBO07_Y_SURFACE * (z - z_if) / VBO_L_PMMA

# The concentrations [PMMA, light solid] of a fully dense mixture at mass
# fraction `y` of the light solid.
function _vbo07_concentrations(y)
    total = 1 / ((1 - y) / VBO_RHO_PMMA + y / VBO07_RHO_LIGHT)
    return (1 - y) * total, y * total
end

function _vbo07_material(;
    substrate_alpha = 8.0e6 / VBO_RHO_KW,
    pmma_alpha = 2870.0 / VBO_RHO_PMMA,
)
    rx = Reaction(:pmma_gasification, 1 => 4; A = 1.35e11, E = 164e3, h = 8.7e5, n = 1.0)
    rx_light =
        Reaction(:light_gasification, 2 => 4; A = 1.35e10, E = 164e3, h = 8.7e5, n = 1.0)
    return Material(
        name = :PMMA_mixed,
        components = (
            _vbo_pmma(; pmma_alpha),
            _vbo_pmma(; name = :PMMA_LIGHT, ρ = VBO07_RHO_LIGHT, pmma_alpha),
            _vbo_kaowool(substrate_alpha),
            _vbo_gas(),
        ),
        reactions = (rx, rx_light),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
end

function _vbo07_problem(n_pmma::Int; t_end = 500.0)
    dz = VBO_L_PMMA / n_pmma
    n_kw = round(Int, VBO_L_KW / dz)
    n = n_pmma + n_kw
    mesh = core_mesh(
        create_uniform_mesh,
        n * dz,
        n,
        Val(4);
        T_initial = VBO_T0,
        ξ_initial = (VBO_RHO_PMMA, 0.0, 0.0, 0.0),
    )
    z_if = n_kw * dz
    ξ0 = [
        z -> z > z_if ? _vbo07_concentrations(_vbo07_fraction(z, z_if))[1] : 0.0,
        z -> z > z_if ? _vbo07_concentrations(_vbo07_fraction(z, z_if))[2] : 0.0,
        z -> z > z_if ? 0.0 : VBO_RHO_KW,
        z -> 0.0,
    ]
    return _vbo_deck(mesh, _vbo07_material(), ξ0; t_end)
end

function _verification_case_vbo07()
    id = "V-BO-07"
    n_default = 60
    function qois(sol)
        l_kw = first(sol.extras.thickness_history) - VBO_L_PMMA
        m = _vbo_measures(sol; l_kw, id)
        # The whole layer leaves: the run ends at the Kaowool's mass.
        exact = 100 * (1 - VBO_RHO_KW * l_kw / first(sol.extras.total_mass))
        out = _vbo_common_qois(
            m;
            id,
            peak_tol = 0.03,
            peak_ref = VBO07_PEAK_REF,
            mass_loss_exact = exact,
        )
        append!(
            out,
            _vbo_timing_qois(m; t_peak_ref = VBO07_T_PEAK_REF, t_bo_ref = VBO07_T_BO_REF),
        )
        VERIFICATION_CORE == "new" &&
            push!(out, _vbo_invocation_qoi(sol, id; merged_limit = nothing))
        return out
    end
    return VerificationCase(
        id = id,
        title = "Burnout of a heterogeneous PMMA layer of mixed composition",
        tags = [:ale, :burnout, :depletion, :beer_lambert, :coupled],
        references = [
            "MaCFP-2 PMMA gasification case (NIST), 50 kW/m², with a second solid mixed in",
        ],
        build = n -> _vbo07_problem(n),
        solve_fn = _vbo_solve_fn(id),
        solve_kwargs = VBO_SOLVE_KWARGS,
        qois = qois,
        ci = CIConfig(n_cells = n_default, budget_s = 20.0),
        report = (curves = sol -> _vbo_history_curves(id, sol),),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo07)
