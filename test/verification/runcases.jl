"""
Verification-case harness (assert mode).

Runs every case in `test/verification/cases/` against its exact solution at
the CI resolution, with pinned tolerances and a convergence-order spot check.
The shared case definitions live in `casedefs.jl` (Test-free), which the
report-mode runner (`scripts/verification/run_reports.jl`) also consumes —
assert mode deliberately does no plotting and no artifact I/O.

Discovery: `test/runtests.jl` auto-includes `.jl` files directly under
`test/verification/` (non-recursive), so this file is the single entry point;
`casedefs.jl`, `cases/`, and `exact/` are only reached through the includes
below.
"""

using Test

include(joinpath(@__DIR__, "casedefs.jl"))

if !isdefined(@__MODULE__, :run_case_ci)
    # Set PYROLYSIS_VERIFICATION_PRINT_QOIS=1 to print every observed QoI and
    # convergence rate — the workflow for (re)pinning tolerances at 3× observed.
    _print_qois() = get(ENV, "PYROLYSIS_VERIFICATION_PRINT_QOIS", "0") == "1"

    function run_case_ci(case::VerificationCase)
        @testset "$(case.id) — $(case.title)" begin
            wall = @elapsed begin
                problem = case.build(case.ci.n_cells)
                sol = _case_solve(case, problem; case.solve_kwargs...)
                @test sol.retcode == :Success
                for q in case.qois(sol)
                    @testset "$(q.name)" begin
                        _print_qois() && println(
                            "[QoI] $(case.id) | $(q.name) | value=$(q.value) " *
                            "exact=$(q.exact) error=$(q.error) tol=$(q.tol)",
                        )
                        if !(q.error < q.tol)
                            @info "QoI out of tolerance" case.id q.name q.value q.exact q.error q.tol q.comment
                        end
                        @test q.error < q.tol
                    end
                end
                conv = case.ci.convergence
                if conv !== nothing
                    @testset "spatial convergence (expected order $(conv.expected_order))" begin
                        build = conv.build === nothing ? case.build : conv.build
                        errors = Float64[]
                        hs = Float64[]
                        for n in conv.n_cells
                            s = _case_solve(case, build(n); conv.solve_kwargs...)
                            push!(errors, conv.error_fn(s))
                            push!(hs, 1.0 / n)
                        end
                        rates = observed_orders(errors, hs)
                        _print_qois() && println(
                            "[QoI] $(case.id) | convergence | errors=$(errors) rates=$(rates)",
                        )
                        if !all(r -> r > conv.min_order, rates)
                            @info "convergence rates" case.id errors rates
                        end
                        @test all(r -> r > conv.min_order, rates)
                        @test mean(rates) > conv.min_mean_order
                    end
                end
            end
            @test wall < 5 * case.ci.budget_s
        end
    end
end

@testset "Verification suite (assert mode)" begin
    include(joinpath(@__DIR__, "exact", "selftests.jl"))

    # invokelatest: the case closures are defined during load (world-age).
    for case in load_verification_cases()
        Base.invokelatest(run_case_ci, case)
    end
end
