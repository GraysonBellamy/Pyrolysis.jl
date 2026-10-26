"""
Every precompile workload runs to its end without a warning. During
precompilation a failing workload only warns, and the configuration it covers
then compiles on first use, so a broken workload would otherwise go unnoticed.
"""

using Test
using Logging
using Pyrolysis

@testset "precompile workloads" begin
    workloads = Pyrolysis._precompile_workloads()
    @test length(workloads) == length(unique(first.(workloads)))
    for (name, run) in workloads
        @testset "$name" begin
            sol = @test_logs min_level = Logging.Warn run()
            if occursin("deactivation", name)
                @test sol.retcode in (:Success, :Terminated)
                @test !isempty(sol.events)
            else
                @test sol.retcode === :Success
            end
        end
    end
end
