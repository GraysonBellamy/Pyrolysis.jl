# =============================================================================
# V-BO-10 — heated burnout without a residue at orders other than one
# =============================================================================
#
# One solid gasifies with no residue under 40 kW/m² on a 10 mm slab, at orders 0,
# 0.5 and 1.2, on Lagrangian cells. Nothing holds a spent cell live, so the column
# burns through to the thickness floor, cell by cell, through the event path; at
# order 1.2, where the cells decay algebraically, at a threshold of 1e-2, since at
# the default none reaches it before the floor. At order 0.5 the threshold's effect
# and the mesh's are bounded: thresholds 1e-2 and 1e-3 agree over the burn, and so
# do 40 and 80 cells, over 5 to 115 s, before the terminal spike of the last cells.
# The two energy forms part there too, by the capacity floor, which is not the
# same relation in each: a bound on that keeps the difference in sight.

const VBO10_L, VBO10_FLUX = 10e-3, 4e4
const VBO10_T_END = 400.0

# The saved interval of the pairs, over the window they compare. Each event is a
# step in the rate, so their peak differences read the time between two runs'
# steps: at the case's 0.5 s they still move with the interval, at this one they
# have settled.
const VBO10_PAIR_DT = 0.025

function _vbo10_problem(N, n; dt_out = 0.5)
    s = SolidComponent(:a, ρ = 500.0, φ = 0.3, c = 1500.0, k = 0.2, ε = 0.9, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1500.0, k = 0.02)
    material = Material(
        name = :VBO10,
        components = (s, gas),
        reactions = (Reaction(:burn, :a => :gas; A = 1e10, E = 1.5e5, h = 3e5, n),),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(VBO10_L, N),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(VBO10_FLUX), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 500.0, z -> 0.0],
        tspan = (0.0, VBO10_T_END),
        output_times = sort!(unique!(vcat(0.0:dt_out:115.0, 0.0:0.5:VBO10_T_END))),
    )
end

const VBO10_SOLVE_KWARGS = (;
    abstol = 1e-9,
    reltol = 1e-7,
    use_ale = true,
    min_thickness = 1e-5,
    dtmin = 1e-14,
    deactivation = true,
)

# Over `a` to `b` s: the largest difference of the gas generation of `test` and
# `ref` over the peak of `ref`, and their integrated difference over the mass the
# reference lost.
function _vbo10_difference(ref, test, a, b)
    window = findall(t -> a <= t <= b, ref.t)
    r = ref.extras.gas_generation_rate[window]
    s = test.extras.gas_generation_rate[window]
    t = ref.t[window]
    d = abs.(s .- r)
    integral = sum((d[j] + d[j+1]) / 2 * (t[j+1] - t[j]) for j = 1:(length(t)-1))
    lost = ref.extras.total_mass[1] - ref.extras.total_mass[window[end]]
    return maximum(d) / maximum(r), integral / lost
end

function _verification_case_vbo10(suffix, n, events; keywords = (;))
    id = "V-BO-10$(suffix)"
    function qois(sol)
        count = length(sol.events)
        ended = sol.outcome === :min_thickness
        out = [
            QoIResult(
                name = "outcome :min_thickness",
                value = ended ? 1.0 : 0.0,
                exact = 1.0,
                error = ended ? 0.0 : 1.0,
                tol = 0.5,
                comment = "the solve ended with outcome :$(sol.outcome)",
            ),
            QoIResult(
                name = "cells through the event path",
                value = Float64(count),
                exact = Float64(events),
                error = abs(count - events),
                tol = 2.5,
                comment = "of 40; the rest hold under 1e-3 of the mass at the floor",
            ),
        ]
        n == 0.5 || return out
        options = VBO10_SOLVE_KWARGS
        pair(N, ε) = suite_solve(
            _vbo10_problem(N, n; dt_out = VBO10_PAIR_DT);
            options...,
            threshold = ε,
        )
        coarse, fine, refined = pair(40, 1e-2), pair(40, 1e-3), pair(80, 1e-3)
        peak, integral = _vbo10_difference(fine, coarse, 5.0, 115.0)
        mesh_peak, mesh_integral = _vbo10_difference(refined, fine, 5.0, 115.0)
        form(energy_state) = solve(
            _vbo10_problem(40, n; dt_out = VBO10_PAIR_DT);
            options...,
            threshold = 1e-3,
            energy_state,
            progress = false,
        )
        form_peak, form_integral =
            _vbo10_difference(form(:temperature), form(:enthalpy), 5.0, 115.0)
        append!(
            out,
            [
                QoIResult(
                    name = "thresholds 1e-2 and 1e-3, peak gas generation",
                    value = peak,
                    error = peak,
                    tol = 0.05,
                    comment = "over 5 to 115 s",
                ),
                QoIResult(
                    name = "thresholds 1e-2 and 1e-3, integrated gas generation",
                    value = integral,
                    error = integral,
                    tol = 0.01,
                    comment = "over the mass the reference lost",
                ),
                QoIResult(
                    name = "40 and 80 cells, peak gas generation",
                    value = mesh_peak,
                    error = mesh_peak,
                    tol = 0.05,
                    comment = "over 5 to 115 s",
                ),
                QoIResult(
                    name = "40 and 80 cells, integrated gas generation",
                    value = mesh_integral,
                    error = mesh_integral,
                    tol = 0.02,
                    comment = "over the mass the reference lost",
                ),
                QoIResult(
                    name = "temperature and enthalpy forms, peak gas generation",
                    value = form_peak,
                    error = form_peak,
                    tol = 0.02,
                    comment = "over 5 to 115 s, the capacity floor's difference",
                ),
                QoIResult(
                    name = "temperature and enthalpy forms, integrated gas generation",
                    value = form_integral,
                    error = form_integral,
                    tol = 0.005,
                    comment = "over the mass the reference lost",
                ),
            ],
        )
        return out
    end
    return VerificationCase(
        id = id,
        title = "Heated burnout without a residue, order $n",
        tags = [:ale, :burnout, :depletion, :coupled],
        references = ["mesh and threshold pairs of the same deck"],
        build = N -> _vbo10_problem(N, n),
        solve_kwargs = (; VBO10_SOLVE_KWARGS..., keywords...),
        expected_retcode = :Terminated,
        qois = qois,
        ci = CIConfig(n_cells = 40, budget_s = n == 0.5 ? 120.0 : 40.0),
    )
end

push!(VERIFICATION_CASE_BUILDERS, () -> _verification_case_vbo10("a", 0.0, 39))
push!(VERIFICATION_CASE_BUILDERS, () -> _verification_case_vbo10("b", 0.5, 38))
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _verification_case_vbo10("c", 1.2, 37; keywords = (; threshold = 1e-2)),
)
