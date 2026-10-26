# Native adaptive-time sensitivity on each comparison's interval-count mesh.
# This is accuracy qualification evidence, not a performance measurement.
using TOML
using Statistics: median
using OrdinaryDiffEq: DiscreteCallback
include(joinpath(@__DIR__, "..", "..", "test", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "report_harness.jl"))
const STUDY_CASES = Dict(c.id => c for c in load_verification_cases())

function study_comparison_times(outdir; wanted = "all")
    targets = TOML.parsefile(joinpath(REFERENCE_DIR, "comparison_targets.toml"))["targets"]
    cases = STUDY_CASES
    work, changes = [], []
    mkpath(outdir)
    for target in targets
        id, n, variant = target["case"], target["n_cells"], target["variant"]
        wanted == "all" || id in split(wanted, ',') || continue
        case = cases[id]
        skw = report_solve_kwargs(case)
        solutions = []
        cap = Inf
        for level = 0:2
            steps = Float64[]
            observer = DiscreteCallback(
                (u, t, integrator) -> t > integrator.tprev,
                integrator -> push!(steps, integrator.t - integrator.tprev);
                save_positions = (false, false),
            )
            problem =
                isempty(variant) ? case.build(n) :
                _vcpl03_problem(Symbol(variant); n_cells = n, n_outputs = 33).problem
            kwargs = merge(
                skw,
                (
                    reltol = get(skw, :reltol, 1e-8) / 10.0^level,
                    abstol = get(skw, :abstol, 1e-10) / 10.0^level,
                    dtmax = level == 0 ? get(skw, :dtmax, Inf) : cap / 2.0^level,
                    callback = observer,
                ),
            )
            if haskey(skw, :pressure_tolerance)
                kwargs = merge(
                    kwargs,
                    (pressure_tolerance = skw.pressure_tolerance / 10.0^level,),
                )
            end
            sol = _case_solve(case, problem; kwargs...)
            sol.retcode == case.expected_retcode ||
                error("$id/$variant n=$n level=$level: $(sol.retcode)")
            isempty(steps) && error("$id: no observed accepted steps")
            length(steps) == sol.statistics.accepted ||
                error("$id: step observer/count mismatch")
            level == 0 && (cap = maximum(steps))
            push!(solutions, sol)
            s = sol.statistics
            push!(
                work,
                (
                    id,
                    variant,
                    n,
                    level,
                    kwargs.reltol,
                    kwargs.abstol,
                    kwargs.dtmax,
                    s.accepted,
                    s.rejected,
                    minimum(steps),
                    median(steps),
                    maximum(steps),
                    sum(steps)/length(steps),
                    s.rhs_evaluations,
                    s.jacobians,
                    s.linear_solves,
                ),
            )
            @info "adaptive time study" id variant n level accepted=s.accepted rejected=s.rejected
        end
        for level = 1:2
            a, b = solutions[level], solutions[level+1]
            a.t == b.t || error("$id: saved observation times changed")
            observables = [
                ("temperature_K", a.T, b.T),
                ("bulk_concentrations_kg_m3", a.ξ, b.ξ),
                ("surface_temperature_K", a.surface_T, b.surface_T),
                ("thickness_m", a.extras.thickness_history, b.extras.thickness_history),
            ]
            if !isempty(a.gas.species)
                push!(observables, ("pressure_Pa", a.gas.pressure, b.gas.pressure))
                push!(
                    observables,
                    ("gas_inventory_kg", a.gas.inventories, b.gas.inventories),
                )
            end
            for (label, aa, bb) in observables
                all(isfinite, aa) && all(isfinite, bb) || error("$id: nonfinite $label")
                push!(
                    changes,
                    (
                        id,
                        variant,
                        n,
                        level-1,
                        level,
                        label,
                        maximum(abs.(aa .- bb)),
                        length(aa),
                    ),
                )
            end
        end
        # Keep completed results if a later case fails.
        write_csv(
            joinpath(outdir, "comparison_julia_steps.csv"),
            [
                "case",
                "variant",
                "cells",
                "level",
                "reltol",
                "abstol",
                "dtmax_s",
                "accepted",
                "rejected",
                "dt_min_s",
                "dt_median_s",
                "dt_max_s",
                "dt_mean_s",
                "rhs_evaluations",
                "jacobians",
                "linear_solves",
            ],
            work,
        )
        write_csv(
            joinpath(outdir, "comparison_julia_time_changes.csv"),
            [
                "case",
                "variant",
                "cells",
                "from_level",
                "to_level",
                "observable",
                "max_change",
                "samples",
            ],
            changes,
        )
    end
end

outdir =
    length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "..", "docs", "guide", "verification", "assets")
Base.invokelatest(
    study_comparison_times,
    outdir;
    wanted = length(ARGS) >= 2 ? ARGS[2] : "all",
)
