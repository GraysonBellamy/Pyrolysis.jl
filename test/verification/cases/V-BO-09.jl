# =============================================================================
# V-BO-09 — burnout through the event path on the conversion
# =============================================================================
#
# A column at a rate that reads no temperature, with reactant A on `basis = 1`
# at the order of the case, its initial conversion graded from 0.99 at the back
# to 0.05 at the exposed face, and a filler B stripped at first order. Every
# material point has a closed form: x = [f^(1−n) − (1−n) k t]₊^(1/(1−n)) and the
# filler (1 − f) exp(−k_B t), so the column's condensed mass and its gas
# generation are known. The burn takes the cells one by one from the exposed
# face, each below the threshold on its reactive inventory, each through a merge
# into the live cell beneath it.
#
# The column must end on `min_thickness` or `burned_out`, with the condensed mass
# on the exact curve and the gas generation on the exact one after the filler is
# gone, its error the transfer's `ε^n` of one cell's rate over the cells burning.

const VBO09_L, VBO09_T, VBO09_ρ = 10e-3, 600.0, 500.0
const VBO09_K_B = 2.0
const VBO09_F_LOW, VBO09_F_HIGH = 0.05, 0.99
const VBO09_THRESHOLD = 1e-4

function _vbo09_problem(n, N, k, t_end)
    A = SolidComponent(:A, ρ = VBO09_ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    B = SolidComponent(:B, ρ = VBO09_ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :VBO09,
        components = (A, B, gas),
        reactions = (
            Reaction(:burn, 1 => 3; A = k, E = 0.0, h = 0.0, n = n, basis = 1.0),
            Reaction(:strip, 2 => 3; A = VBO09_K_B, E = 0.0, h = 0.0),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1e-3, conversion_threshold = 1e-6),
    )
    f(z) = VBO09_F_HIGH - (VBO09_F_HIGH - VBO09_F_LOW) * z / VBO09_L
    return PyrolysisProblem(
        mesh = create_uniform_mesh(VBO09_L, N),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = VBO09_T)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> VBO09_T,
        ξ_initial = [z -> VBO09_ρ * f(z), z -> VBO09_ρ * (1 - f(z)), z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:(t_end/1000):t_end),
    )
end

# The exact mass and bare rate per unit initial condensed mass, by the midpoint
# rule over the column's initial composition.
function _vbo09_exact(n, k, t; points = 40_000)
    M, r = 0.0, 0.0
    for q = 1:points
        f = VBO09_F_LOW + (VBO09_F_HIGH - VBO09_F_LOW) * (q - 0.5) / points
        p = f^(1 - n) - (1 - n) * k * t
        x = p > 0 ? p^(1 / (1 - n)) : 0.0
        b = (1 - f) * exp(-VBO09_K_B * t)
        M += x + b
        r += (x > 0 ? k * x^n : 0.0) + VBO09_K_B * b
    end
    return M / points, r / points
end

function _verification_case_vbo09(suffix, n, k, t_end, mass_tol, rate_tol, events)
    id = "V-BO-09$(suffix)"
    function qois(sol)
        condensed = sol.m[1, :, :] .+ sol.m[2, :, :]
        mass = vec(sum(condensed, dims = 1)) ./ sum(condensed[:, 1])
        rate = sol.extras.gas_generation_rate ./ (VBO09_ρ * VBO09_L)
        pairs = [_vbo09_exact(n, k, t) for t in sol.t]
        exact_M, exact_r = first.(pairs), last.(pairs)
        after = sol.t .>= 6.0
        mass_error = maximum(abs.(mass .- exact_M))
        rate_error = maximum(abs.(rate[after] .- exact_r[after])) / maximum(exact_r[after])
        ended = sol.outcome in (:min_thickness, :burned_out)
        return [
            QoIResult(
                name = "condensed mass over its initial value, largest error",
                value = mass_error,
                error = mass_error,
                tol = mass_tol,
                comment = "three times the observed",
            ),
            QoIResult(
                name = "gas generation after the filler, largest error over the exact peak",
                value = rate_error,
                error = rate_error,
                tol = rate_tol,
                comment = "three times the observed",
            ),
            QoIResult(
                name = "cells through the event path",
                value = Float64(length(sol.events)),
                exact = Float64(events),
                error = abs(length(sol.events) - events),
                tol = 3.5,
                comment = "the column ends on the thickness with a few cells left",
            ),
            QoIResult(
                name = "outcome :min_thickness or :burned_out",
                value = ended ? 1.0 : 0.0,
                exact = 1.0,
                error = ended ? 0.0 : 1.0,
                tol = 0.5,
                comment = "the solve ended with outcome :$(sol.outcome)",
            ),
        ]
    end
    return VerificationCase(
        id = id,
        title = "Burnout through the event path on the conversion, order $n",
        tags = [:ale, :burnout, :depletion, :kinetics],
        references = ["closed-form decay of each material point"],
        build = N -> _vbo09_problem(n, N, k, t_end),
        solve_kwargs = (
            use_ale = true,
            deactivation = true,
            threshold = VBO09_THRESHOLD,
            reltol = 1e-6,
            abstol = 1e-8,
        ),
        expected_retcode = :Terminated,
        qois = qois,
        ci = CIConfig(n_cells = 40, budget_s = 40.0),
    )
end

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _verification_case_vbo09("a", 0.5, 0.02, 110.0, 3 * 2.0e-5, 3 * 4.4e-4, 36),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _verification_case_vbo09("b", 1.2, 0.2, 300.0, 3 * 5.2e-6, 3 * 1.0e-5, 19),
)
