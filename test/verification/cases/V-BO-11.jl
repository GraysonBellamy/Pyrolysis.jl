# =============================================================================
# V-BO-11 — a cell that dies under a live one at order one
# =============================================================================
#
# The PMMA cone of examples/06 on Lagrangian cells with deactivation. The tail a
# merge hands its receiver delays it past its inward neighbour, so a cell crosses
# the threshold under a live cell: at 476 s on this deck. It goes into the cell
# beneath it, the live cells above move down, and the solve runs on. The mass
# loss rate stays within 0.5 % of the peak of the corrected default's, which
# differs only by the correction's redistribution of material.

const VBO11_SECONDS = 600.0

# examples/06 at `N` cells over the case's window.
function _vbo11_problem(N)
    lines = readlines(
        joinpath(@__DIR__, "..", "..", "..", "examples", "06_cone_calorimeter.jl"),
    )
    a = findfirst(startswith("pmma_cp = CallableProperty("), lines)
    b = findfirst(startswith("solution = solve("), lines)
    out = String[]
    skip = false
    for line in lines[a:(b-1)]
        if skip
            startswith(line, ")") && (skip = false)
            continue
        end
        if startswith(line, "println(")
            skip = !endswith(rstrip(line), ")")
            continue
        end
        startswith(line, "flush(") && continue
        push!(out, line)
    end
    m = Module(:VBO11Deck)
    Core.eval(m, :(using Pyrolysis))
    Core.eval(m, :(const SMOKE = false))
    Core.eval(m, Meta.parseall(join(out, "\n")))
    return Core.eval(
        m,
        :(PyrolysisProblem(
            mesh = uniform_mesh(5.8e-3, $N),
            material = problem.material,
            bc_set = problem.bc_set,
            T_initial = problem.T_initial,
            ξ_initial = problem.ξ_initial,
            tspan = (0.0, $(VBO11_SECONDS)),
            dt_initial = problem.dt_initial,
            output_times = collect(0.0:1.0:($(VBO11_SECONDS))),
        )),
    )
end

const VBO11_SOLVE_KWARGS = (;
    abstol = 1e-11,
    reltol = 1e-9,
    pressure_tolerance = 1e-4,
    dtmin = 1e-14,
    use_ale = true,
    min_thickness = 0.01e-3,
    deactivation = true,
)

function _verification_case_vbo11()
    id = "V-BO-11"
    function qois(sol)
        problem = Base.invokelatest(_vbo11_problem, 72)
        ref = Base.invokelatest(suite_solve, problem; VBO11_SOLVE_KWARGS...)
        n = min(length(sol.t), length(ref.t))
        sol.t[1:n] == ref.t[1:n] || error("$id: the two solves saved different times")
        window = findall(>=(30.0), sol.t[1:n])
        r = ref.extras.gas_generation_rate[window]
        s = sol.extras.gas_generation_rate[window]
        peak = maximum(abs.(s .- r)) / maximum(r)
        interior = count(e -> e.cell != 72 - (e.invocation - 1), sol.events)
        return [
            QoIResult(
                name = "gas generation against the corrected default, largest over the peak",
                value = peak,
                error = peak,
                tol = 0.005,
                comment = "from 30 s to the end",
            ),
            QoIResult(
                name = "cells gone from under a live cell",
                value = Float64(interior),
                error = interior >= 1 ? 0.0 : 1.0,
                tol = 0.5,
                comment = "the first at 476 s on this deck",
            ),
            QoIResult(
                name = "the solve reaches its end",
                value = last(sol.t),
                exact = VBO11_SECONDS,
                error = abs(last(sol.t) - VBO11_SECONDS),
                tol = 1e-6,
            ),
        ]
    end
    return VerificationCase(
        id = id,
        title = "A cell that dies under a live one at order one",
        tags = [:ale, :burnout, :depletion, :coupled],
        references = ["the corrected default of the same deck"],
        build = N -> Base.invokelatest(_vbo11_problem, N),
        # The deck's closures are defined when the case builds it.
        solve_fn = (problem; kwargs...) ->
            Base.invokelatest(suite_solve, problem; kwargs...),
        solve_kwargs = (;
            VBO11_SOLVE_KWARGS...,
            correction = Pyrolysis.Closures.NoCorrection(),
        ),
        qois = qois,
        ci = CIConfig(n_cells = 72, budget_s = 120.0),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo11)
