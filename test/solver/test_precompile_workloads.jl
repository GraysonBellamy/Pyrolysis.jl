"""
Every precompile workload runs to its end. During precompilation a failing
workload only warns, and the configuration it covers then compiles on first
use, so a broken workload would otherwise go unnoticed.
"""

using Test
using Pyrolysis

@testset "precompile workloads" begin
    workloads = Pyrolysis._precompile_workloads()
    @test length(workloads) == length(unique(first.(workloads)))
    for (name, run) in workloads
        @testset "$name" begin
            sol = run()
            @test sol.retcode in (:Success, :Terminated)
            if occursin("depletion", name)
                @test !isempty(sol.extras.park_times)
            end
        end
    end
end
