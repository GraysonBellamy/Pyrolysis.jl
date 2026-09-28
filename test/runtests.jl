using Test
using Pyrolysis

# Each file under jacobian/ is evaluated in its own throwaway module: the
# operator/structured tests define top-level structs and consts that would
# clash if dumped into one shared namespace. The module name is derived from
# the path so new files are picked up automatically.
function _include_isolated(relpath::AbstractString)
    name = Symbol(replace(relpath, r"[/\\.\-]" => "_"))
    mod = Module(name)
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Core.eval(mod, :(eval(x) = Core.eval($mod, x)))
    Base.include(mod, joinpath(@__DIR__, relpath))
    return mod
end

function _test_project_for_subprocess()
    project = Base.active_project()
    project_dir = dirname(project)
    if occursin("BenchmarkTools", read(project, String))
        return project_dir
    end
    return @__DIR__
end

function _run_isolated_test_script(relpath::AbstractString)
    project_dir = _test_project_for_subprocess()
    script = joinpath(@__DIR__, relpath)
    cmd = `$(Base.julia_cmd()) --project=$project_dir --startup-file=no --check-bounds=yes $script`
    return success(cmd)
end

_jl_files(dir) = sort(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, dir))))

@testset "Pyrolysis.jl" begin
    # Domain suites — shared namespace; shared fixtures come from helpers.jl
    # (idempotent under repeated include via its own load guard).
    for dir in [
        "boundary_conditions",
        "discretization",
        "solver",
        "verification",
        "ale",
        "adaptivity",
        "radiation",
        "geometry",
        "materials",
        "physics",
        "residual",
    ]
        for f in _jl_files(dir)
            include(joinpath(@__DIR__, dir, f))
        end
    end

    # Jacobian suites — one module each.
    for f in _jl_files("jacobian")
        _include_isolated(joinpath("jacobian", f))
    end

    # Performance invariants are process-state sensitive, especially the
    # allocation probes. Run the perf driver in a fresh Julia process so
    # `Pkg.test()` exercises the same hermetic path as the CI perf job.
    @testset "performance invariants (isolated)" begin
        @test _run_isolated_test_script("perf/runtests.jl")
    end
end
