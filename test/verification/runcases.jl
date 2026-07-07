"""
Verification-case harness (assert mode).

Runs every case in `test/verification/cases/` against its exact solution at
the CI resolution, with pinned tolerances and a convergence-order spot check.
The same `VerificationCase` metadata is designed to be consumed later by the
report-mode runner (`scripts/verification/`, Phase 2) and the guide chapter
generator — assert mode deliberately does no plotting and no artifact I/O.

Coordinate convention (matches `test_analytical.jl`): z = 0 bottom/substrate,
z = L top/exposed surface; exact solutions use depth x = L − z.

Discovery: `test/runtests.jl` auto-includes `.jl` files directly under
`test/verification/` (non-recursive), so this file is the single entry point;
`cases/` and `exact/` are only reached through the includes below.
"""

using Test
using Pyrolysis
using Pyrolysis.Solver: get_cell_positions
using Statistics: mean
using SpecialFunctions: erfc, erfcx

if !isdefined(@__MODULE__, :VerificationCase)
    include(joinpath(@__DIR__, "exact", "exact_solutions.jl"))

    """
    One quantity of interest evaluated from a solved case: `error < tol` is the
    assertion. `comment` records tolerance provenance (observed error at pin
    time, per cross-cutting requirement 5 of the suite plan).
    """
    Base.@kwdef struct QoIResult
        name::String
        value::Float64
        exact::Float64 = NaN
        error::Float64
        tol::Float64
        comment::String = ""
    end

    """
    Spatial convergence spot check: solve `build(n)` for each `n ∈ n_cells`,
    reduce each solution to a scalar error norm via `error_fn`, and assert the
    observed order between consecutive resolutions. `build === nothing` reuses
    the case's own builder. Δz ∝ 1/n is used for the order fit, so `n_cells`
    must refine the same mesh family.
    """
    Base.@kwdef struct ConvergenceSpec
        n_cells::Vector{Int}
        error_fn::Function
        expected_order::Float64
        min_order::Float64
        min_mean_order::Float64
        solve_kwargs::NamedTuple = (;)
        build::Union{Nothing,Function} = nothing
    end

    Base.@kwdef struct CIConfig
        n_cells::Int
        budget_s::Float64            # recorded target; asserted at 5× for CI variance
        convergence::Union{Nothing,ConvergenceSpec} = nothing
    end

    """
    Single-source verification-case definition (assert mode consumes `ci`;
    report mode will consume `report`). `build(n_cells) -> PyrolysisProblem`,
    `qois(sol) -> Vector{QoIResult}`.
    """
    Base.@kwdef struct VerificationCase
        id::String
        title::String
        tags::Vector{Symbol}
        references::Vector{String}
        build::Function
        solve_kwargs::NamedTuple = (;)
        qois::Function
        ci::CIConfig
        report::NamedTuple = (;)
    end

    observed_orders(errors, hs) =
        [log(errors[i-1] / errors[i]) / log(hs[i-1] / hs[i]) for i = 2:length(errors)]

    const VERIFICATION_CASE_BUILDERS = Function[]

    function run_case_ci(case::VerificationCase)
        @testset "$(case.id) — $(case.title)" begin
            wall = @elapsed begin
                problem = case.build(case.ci.n_cells)
                sol = solve(problem; progress = false, case.solve_kwargs...)
                @test sol.retcode == :Success
                for q in case.qois(sol)
                    @testset "$(q.name)" begin
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
                            s = solve(build(n); progress = false, conv.solve_kwargs...)
                            push!(errors, conv.error_fn(s))
                            push!(hs, 1.0 / n)
                        end
                        rates = observed_orders(errors, hs)
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

    empty!(VERIFICATION_CASE_BUILDERS)
    for f in sort(filter(endswith(".jl"), readdir(joinpath(@__DIR__, "cases"))))
        include(joinpath(@__DIR__, "cases", f))
    end
    for case_builder in VERIFICATION_CASE_BUILDERS
        run_case_ci(case_builder())
    end
end
